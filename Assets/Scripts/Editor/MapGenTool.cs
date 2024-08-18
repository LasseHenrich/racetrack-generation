using UnityEngine;
using UnityEditor;
using UnityEngine.SceneManagement;
using UnityEditor.SceneManagement;
using System.Collections.Generic;
using System;
using System.Linq;
using UnityEngine.UIElements;
using System.Runtime.Serialization.Formatters.Binary;
using System.IO;

public class MapGenTool : ToolWindow
{
    protected override string ToolLabel => "Map Generation";

    Scene scene;
    GameObject obj; // As we can't serialize, we store everything in an object in the scene

    #region Curve

    #region Polyline

    #region Exposed

    public GenMode genMode = GenMode.Circular;
    public RepulsionType repulsionType = RepulsionType.Sobolev;

    bool curveClosed = true;
    bool deacObsAfterScaling = true;
    bool rotateAfterScaling = false;
    bool noRepulsionAfterScaling = false;
    float lsStepThreshold = 1e-15f;
    float energyThreshold = 0f;
    bool showBezier = true;
    bool showBezierHandles = true;
    bool showPolyLine = true;
    bool showPolyPoints = true;
    bool showObstacles = true;
    bool useBarnesHut = true;
    bool useBackproj = true;
    bool runningLineSearch = false;

    [SerializeField]
    GenModeConfig genModeConfig;

    public EnergyCurve_EditorConfig Config
    {
        get
        {
            return new EnergyCurve_EditorConfig(curveClosed, lengthScale, deacObsAfterScaling, rotateAfterScaling, noRepulsionAfterScaling, lsStepThreshold, energyThreshold, genModeConfig, ObstacleList, numObstacles, PotentialList, ConstraintList);
        }
    }

    #region Constraints

    float lengthScale = 6;

    [SerializeField]
    Dict_ConstraintType_Bool usingConstraint = new();

    List<ConstraintType> ConstraintList
    {
        get
        {
            return usingConstraint.Where(x => x.Value).Select(x => x.Key).ToList();
        }
    }

    #endregion

    #region Potentials

    Dict_PotentialType_Bool usingPotential = new();
    Dict_PotentialType_PotentialConfig potentialDict = new();
    List<PotentialConfig> PotentialList
    {
        get
        {
            return potentialDict.Where(x => usingPotential.ContainsKey(x.Key) && usingPotential[x.Key] == true).Select(x => x.Value).ToList();
        }
    }

    #endregion

    #region Obstacles

    bool overrideObstacles = false;
    int numObstacles = 10;
    readonly List<ObstacleConfig> obstacleList = new();
    List<ObstacleConfig> ObstacleList
    {
        get
        {
            return overrideObstacles ? obstacleList : null;
        }
    }

    #endregion

    #endregion

    int numStuckIterations;
    const int stepLimit = -1;
    int ls_currentStep;

    int steps;

    int subdivideCount;
    const int subdivideLimit = 2;

    List<CurveVertex> PolyPoints { get { return curve.verts; } }

    EnergyCurve curve
    {
        get
        {
            if (obj == null)
                InitObj();
            return obj.GetComponent<EnergyCurveW>().curve;
        }
        set
        {
            obj.GetComponent<EnergyCurveW>().curve = value;
        }
    }

    #endregion

    #region BezierSpline

    // Spline
    RoadSpline roadSpline;
    float epsilon = 0.2f;
    float psi = 2f;

    // Intersections
    float intersectionPreferredDistance = 0f; //5f

    // Topology
    float _fixedPartLengthMultiplier = 1f;
    float FixedPartLengthMultiplier { get { return _fixedPartLengthMultiplier * 1000; } }
    float crossingExtraSize = 0.2f;

    #endregion

    #endregion

    #region RoadGen

    RoadObjectW _roadObjectWrapper;
    RoadObjectW RoadObjectWrapper
    {
        get
        {
            if (_roadObjectWrapper == null)
                _roadObjectWrapper = obj.GetComponent<RoadObjectW>();
            return _roadObjectWrapper;
        }
    }

    float roadPartLength = 1f;
    float bridgePartLength = 1f;
    float _widthMultiplier = 1f;
    float WidthMultiplier { get { return _widthMultiplier * 1000; } }
    float _heightMultiplier = 1f;
    float HeightMultiplier { get { return _heightMultiplier * 1000; } }
    float crossingShape = 0.5f;

    bool autoUpdateRoad = true;

    string mapObjectName = "MyMap";

    #endregion

    /// <summary>
    /// When opening the window, add our Draw() function to the Unity drawing stack
    /// </summary>
    protected override void OnEnable()
    {
        base.OnEnable();

        SceneView.duringSceneGui += MyUpdate;

        try {
            curve.TryInitVertsEdgesFromSerialized(); 
            curve.S_TryInitObstacles();
        }
        catch(Exception)
        {
            Debug.Log("No curve found, generating...");
            GeneratePolyline();
            Debug.Log("Curve generated");
        }
    }

    /// <summary>
    /// When closing the window, remove our Draw() function from the Unity drawing stack
    /// </summary>
    void OnDisable()
    {
        SceneView.duringSceneGui -= MyUpdate;
    }

    void InitObj()
    {
        Debug.Log("InitObj");
        obj = GameObject.FindGameObjectWithTag("RoadObject");
        if (obj == null)
        {
            obj = new()
            {
                name = "RoadObject",
                tag = "RoadObject"
            };

            obj.AddComponent<EnergyCurveW>();
            obj.AddComponent<RoadObjectW>();
        }
    }

    [MenuItem("Window/Map Gen Tool")]
    public static void ShowWindow()
    {
        GetWindow(typeof(MapGenTool));
    }

    protected override void OnGUI()
    {
        base.OnGUI();

        #region Manual Steps
        CreateSection("Manual Steps", () =>
        {
            CreateLabel("Shape Generation", 16);
            CreateButton("Create Polyline", () => GeneratePolyline());
            CreateButton("Reset Polyline", ResetPolyline);
            CreateButton("Single Step", () => RepulsionUpdate());
            CreateButton("Ten Steps", () => RepulsionUpdate(10));
            CreateCheckbox("Run Shape Generation automatically", ref runningLineSearch);

            GUILayout.Space(5);

            CreateLabel("Spline", 16);
            CreateButton("Generate Bezier Spline", GeneateBezierSpline);
            CreateButton("Add Intersection", () => roadSpline.AddIntersections(intersectionPreferredDistance));

            GUILayout.Space(5);

            CreateLabel("Mesh Generation", 16);
            CreateButton("Scale Spline to fit Width", ScaleCurveToFitWidth);
            CreateButton("Calculate Features along Spline", () => TopologyHandler.GenerateTopologies());
            CreateButton("Generate Road", GenerateRoad);

            GUILayout.Space(5);

            CreateLabel("Saving", 16);
            CreateTextField("File Name", ref mapObjectName);
            CreateButton("Save Map", SaveMap);
            CreateButton("Export Prefab", ExportPrefab);
        });
        #endregion

        GUILayout.Space(10);

        #region Advanced
        CreateSection("Advanced", () =>
        {

            #region Curve Generation
            CreateSubSection("Curve Generation", () =>
            {

                #region GUI
                CreateSubSubSection("GUI", () =>
                {
                    CreateCheckbox("Show Bezier", ref showBezier);
                    CreateCheckbox("Show Bezier Handles", ref showBezierHandles);
                    CreateCheckbox("Show Poly Line", ref showPolyLine);
                    CreateCheckbox("Show Poly Points", ref showPolyPoints);
                    CreateCheckbox("Show Obstacles", ref showObstacles);
                });
                #endregion

                #region Polyline
                CreateSubSubSection("Polyline", () =>
                {

                    genMode = (GenMode)EditorGUILayout.EnumPopup("Gen Mode", genMode);
                    CreateCheckbox("Closed Curve", ref curveClosed);

                    #region GenModeConfig
                    if (genMode == GenMode.Bezier)
                    {
                        if (genModeConfig is not GenModeConfig_Bezier) genModeConfig = new GenModeConfig_Bezier();
                        CreateSlider(ref ((GenModeConfig_Bezier)genModeConfig).numPoints, "Spline Points", 3, 10);
                        CreateSliderFloat(ref ((GenModeConfig_Bezier)genModeConfig).radius, "Radius", 1, 20);
                        CreateSliderFloat(ref ((GenModeConfig_Bezier)genModeConfig).spacing, "Spacing", 0.2f, 2f);
                    }
                    else if (genMode == GenMode.Circular)
                    {
                        if (genModeConfig is not GenModeConfig_Circular) genModeConfig = new GenModeConfig_Circular();
                        CreateSlider(ref ((GenModeConfig_Circular)genModeConfig).numPoints, "Points", 10, 100);
                        CreateSliderFloat(ref ((GenModeConfig_Circular)genModeConfig).radius, "Radius", 1, 20);
                    }
                    #endregion

                    #region Obstacles

                    CreateLabel("");
                    CreateLabel("Obstacles");

                    CreateSlider(ref numObstacles, "Obstacles", 0, 20);
                    CreateCheckbox("Override Obstacles", ref overrideObstacles);

                    if (overrideObstacles)
                    {

                        EditorGUIUtility.labelWidth = 50;
                        EditorGUIUtility.fieldWidth = 20;
                        for (int i = 0; i < numObstacles; i++)
                        {
                            GUILayout.BeginHorizontal();
                            GUILayout.Label("Obstacle " + (i + 1));

                            if (i >= obstacleList.Count) obstacleList.Add(new(weight: 1, radius: 5, numPoints: 20, center: Vector2.zero));

                            ObstacleConfig obs = obstacleList[i];

                            CreateFloatField("Weight", ref obs.weight);
                            CreateIntField("Points", ref obs.numPoints);
                            CreateFloatField("Radius", ref obs.radius);
                            CreateVector2Field("Center", ref obs.center);

                            GUILayout.EndHorizontal();
                        }

                        if (obstacleList.Count > numObstacles)
                        {
                            obstacleList.RemoveRange(Mathf.Max(numObstacles - 1, 0), obstacleList.Count - numObstacles);
                        }

                        // Reset
                        EditorGUIUtility.labelWidth = 0;
                        EditorGUIUtility.fieldWidth = 0;
                    }

                    if (numObstacles > 0)
                    {
                        CreateCheckbox("Deactivate after length scaling", ref deacObsAfterScaling);
                    }

                    #endregion

                    #region Constraints

                    CreateLabel("");
                    CreateLabel("Constraints");
                    foreach (ConstraintType c in Enum.GetValues(typeof(ConstraintType)))
                    {
                        if (!usingConstraint.ContainsKey(c))
                            usingConstraint.Add(c, true);
                        CreateCheckbox_Dict(c.ToString(), usingConstraint, c);
                        if (c == ConstraintType.Length) CreateSliderFloat(ref lengthScale, "Target Scaled Length", 0.5f, 10f);
                    }

                    #endregion

                    #region Potentials

                    CreateLabel("");
                    CreateLabel("Potentials");
                    EditorGUIUtility.labelWidth = 50;
                    foreach (PotentialType pType in Enum.GetValues(typeof(PotentialType)))
                    {
                        GUILayout.BeginHorizontal();

                        if (!usingPotential.ContainsKey(pType))
                            usingPotential.Add(pType, true);
                        CreateCheckbox_Dict(pType.ToString(), usingPotential, pType);

                        if (usingPotential[pType])
                        {
                            if (!potentialDict.ContainsKey(pType))
                                potentialDict.Add(pType, new(pType, 0.1f, true));
                            PotentialConfig pConf = potentialDict.LastOrDefault().Value;
                            CreateFloatField("Weight", ref pConf.weight);
                            CreateCheckbox("Delayed", ref pConf.delayed);
                        }

                        GUILayout.EndHorizontal();
                    }
                    // Reset
                    EditorGUIUtility.labelWidth = 0;

                    #endregion

                    #region Other

                    CreateLabel("");
                    CreateLabel("Other");
                    CreateCheckbox("Rotate after scaling", ref rotateAfterScaling);
                    CreateCheckbox("No repulsion after scaling", ref noRepulsionAfterScaling);
                    CreateFloatField("LS step threshold", ref lsStepThreshold);
                    CreateFloatField("energy threshold", ref energyThreshold);

                    #endregion

                    #region Creation
                    CreateButton("Create Polyline", () => GeneratePolyline());
                    CreateButton("Reset Polyline", ResetPolyline);
                    CreateButton("Print Polyline as 2D Arary", PrintPolyline);
                    #endregion

                    #region Repulsion
                    CreateLabel("");
                    CreateLabel("Repulsion");
                    CreateCheckbox("Use Barnes Hut", ref useBarnesHut);
                    CreateCheckbox("Use Backprojection", ref useBackproj);
                    repulsionType = (RepulsionType)EditorGUILayout.EnumPopup("Repulsion Method", repulsionType);
                    CreateLabel("");
                    CreateLabel("Running");
                    CreateButton("Single Step LS", () => RepulsionUpdate());
                    CreateButton("Ten Steps LS", () => RepulsionUpdate(10));
                    CreateCheckbox("Run LS", ref runningLineSearch);

                    #endregion

                });
                #endregion

                #region Road Spline
                CreateSubSubSection("Bezier Spline", () =>
                {
                    CreateFloatField("Epsilon", ref epsilon);
                    CreateFloatField("Psi", ref psi);
                    CreateButton("Generate Bezier Spline", GeneateBezierSpline);
                    CreateButton("Add Intersection", () => roadSpline.AddIntersections(intersectionPreferredDistance));
                    CreateFloatField("Fixed Part Length Multiplier", ref _fixedPartLengthMultiplier);
                    CreateFloatField("Crossing Extra Size", ref crossingExtraSize);
                    CreateButton("Scale curve to fit width", ScaleCurveToFitWidth);
                    CreateButton("Add Features", () => TopologyHandler.GenerateTopologies());
                });
                #endregion
            });
            #endregion

            #region Road
            CreateSubSection("Road Object", () =>
            {
                // Note that the "name" arguments in the DoFoldout() calls also ensures
                // the hashing in CreateSubSubSection to work correctly

                CreateSubSubSection("Road Config", () =>
                {
                    DoFoldout("Road", RoadObjectWrapper.road);
                });

                CreateSubSubSection("BridgeAscent Config", () =>
                {
                    DoFoldout("BridgeAscent", RoadObjectWrapper.road);
                });

                CreateSubSubSection("Bridge Config", () =>
                {
                    DoFoldout("Bridge", RoadObjectWrapper.road);
                });

                CreateSubSubSection("BridgeDescent Config", () =>
                {
                    DoFoldout("BridgeDescent", RoadObjectWrapper.road);
                });

                CreateSubSubSection("Crossing Config", () =>
                {
                    DoFoldout("Crossing", RoadObjectWrapper.road);
                });

                CreateSubSubSection("Ramp Config", () =>
                {
                    DoFoldout("Ramp", RoadObjectWrapper.road);
                });

                void DoFoldout(string name, MeshTopologyEditorConfig metc)
                {
                    metc.mesh = (Mesh)EditorGUILayout.ObjectField(name + " Mesh", metc.mesh, typeof(Mesh), allowSceneObjects: false);
                    CreateSlider(ref metc.numMaterials, name + " Materials", 0, 5);
                    CreateMaterialSelection(metc.materials, metc.numMaterials);
                }

                void CreateMaterialSelection(List<Material> matList, int numMats)
                {
                    for (int i = 0; i < numMats; i++)
                    {
                        if (i >= matList.Count) matList.Add(new(Shader.Find("Specular")));
                        matList[i] = (Material)EditorGUILayout.ObjectField("Material " + (i + 1), matList[i], typeof(Material), allowSceneObjects: false);
                    }

                    if (matList != null && matList.Count > numMats)
                    {
                        matList.RemoveRange(Mathf.Max(numMats - 1, 0), matList.Count - numMats);
                    }
                }

                CreateFloatField("Width Multiplier", ref _widthMultiplier);
                CreateFloatField("Height Multiplayer", ref _heightMultiplier);
                CreateFloatField("Road Part Length", ref roadPartLength);
                CreateFloatField("Bridge Part Length", ref bridgePartLength);
                CreateSliderFloat(ref crossingShape, "Crossing Shape", 0f, 1f);
                CreateCheckbox("Auto Update Road", ref autoUpdateRoad);
                CreateButton("Generate Road", GenerateRoad);
            });
            #endregion

            #region Operations
            CreateSubSection("Operations", () =>
            {
                //CreateButton("Start Environment", StartEnvironment);
                //CreateButton("Stop Environment", StopEnvironment);
                CreateButton("Load Map", LoadMap);
                CreateTextField("File Name", ref mapObjectName);
                CreateButton("Save Map", SaveMap);
                CreateButton("Export Prefab", ExportPrefab);
            });
            #endregion

        });
        #endregion

        EditorGUILayout.EndScrollView();
        GUILayout.EndVertical();

        UpdateTopologyHandler();
    }

    private void UpdateTopologyHandler()
    {
        TopologyHandler.OnUpdateProperties(
            roadPartLength: roadPartLength,
            roadPartWidth: RoadPartWidth(),
            bridgeAscentLength: BridgeAscentLength(),
            bridgeDescentLength: BridgeDescentLength(),
            crossingExtraSize: crossingExtraSize,
            rampLength: RampLength(),
            roadSpline: roadSpline
        );
    }

    void MyUpdate(SceneView sceneView)
    {
        if (runningLineSearch)
        {
            RepulsionUpdate();
        }

        try
        {
            Draw();
        }
        catch (Exception)
        {
            Debug.Log("No curve found, generating...");
            GeneratePolyline();
            Debug.Log("Curve generated");
        }

        curve.SerializeVertsEdgesPositions();
        curve.S_SerializeObstaclePositions();
    }

    void Draw()
    {
        //if (curve.verts != null) Debug.Log(curve.verts.Count);

        #region Colors

        Color controlHandleColor = new(0.2f, 0.5f, 0.2f);
        Color anchorHandleColor = new(255f / 255f, 225f / 255f, 67f / 255f);
        Color lineColor = new(0, 0, 0);
        Color splineColor = new(0.2f, 0.2f, 0.5f);
        Color polyColor = new(0.4f, 0.6f, 1f);

        float bezierWidth =  8f; // 20f
        float bezierLineWidth = 3f; // 5f
        float bezierHandleWidth = 1.4f;//1.0f;//0.6f;
        float polylineWidth = 8f;//4f;
        float polylineDiscWidth = 0.15f;
        float obstacleWidth = 5f;//2f;
        float obstacleDiscWidth = 0.15f;

        #endregion

        #region Polyline

        #region Polyline
        if (PolyPoints != null)
        {
            Handles.color = polyColor;
            if (showPolyLine)
            {
                List<Vector3> polyPoints3D = new List<Vector3>();
                for (int i = 0; i < PolyPoints.Count; i++) 
                    polyPoints3D.Add(xzToX0Z(PolyPoints[i].Position()));
                if (curveClosed) polyPoints3D.Add(xzToX0Z(PolyPoints[0].Position())); // To draw a closed curve, we need to add the first point again
                Handles.DrawAAPolyLine(EditorGUIUtility.whiteTexture, polylineWidth, polyPoints3D.ToArray());
            }
            if (showPolyPoints)
            {
                for (int i = 0; i < PolyPoints.Count; i++)
                {
                    Handles.DrawSolidDisc(xzToX0Z(PolyPoints[i].Position()), Vector3.up, polylineDiscWidth);
                }
            }
        }
        #endregion

        #region Obstacles
        if (showObstacles && curve != null && curve.obstacles != null)
        {
            Handles.color = Color.red;
            foreach (Obstacle obs in curve.obstacles)
            {
                if (!obs.IsEnabled)
                    continue;

                List<Vector3> obstaclePoints3D = new List<Vector3>();
                foreach (CurveVertex vert in obs.verts)
                    obstaclePoints3D.Add(xzToX0Z(vert.Position()));
                obstaclePoints3D.Add(xzToX0Z(obs.verts[0].Position())); // To draw a closed curve, we need to add the first point again
                if (showPolyLine)
                {
                    Handles.DrawAAPolyLine(EditorGUIUtility.whiteTexture, obstacleWidth, obstaclePoints3D.ToArray());
                }
                if (showPolyPoints)
                {
                    for (int i = 0; i < obstaclePoints3D.Count; i++)
                    {
                        Handles.DrawSolidDisc(obstaclePoints3D[i], Vector3.up, obstacleDiscWidth);
                    }
                }
            }

        }
        #endregion

        #endregion

        #region Bezier Spline
        if (roadSpline != null && showBezier)
        {
            DrawSpline(roadSpline, autoUpdateRoad ? RegenerateRoad : null);
        }
        #endregion

        void DrawSpline(Spline spline, Action OnUpdated = null)
        {
            if (spline == null)
                return;

            #region Lines and Bezier
            Handles.color = lineColor;
            for (int i = 0; i < spline.NumSegments; i++) 
            {
                // Sometimes the spline loses its points, so we need to recalculate them
                try
                {
                    spline.GetPointsInSegment(i);
                }
                catch (Exception)
                {
                    spline.CalculateAllPointsInSegments();
                }
                
                Vector2[] points = spline.GetPointsInSegment(i);

                Vector3[] points3D = new Vector3[points.Length];
                for (int p = 0; p < points.Length; p++)
                {
                    points3D[p] = new Vector3(points[p].x, 0, points[p].y);
                }

                if (showBezierHandles)
                {
                    Handles.DrawLine(points3D[1], points3D[0], bezierLineWidth);
                    Handles.DrawLine(points3D[2], points3D[3], bezierLineWidth);
                }
                Handles.DrawBezier(points3D[0], points3D[3], points3D[1], points3D[2], splineColor, Texture2D.whiteTexture, bezierWidth);
            }
            #endregion

            #region Handles
            if (showBezierHandles)
            {
                for (int i = 0; i < spline.NumPoints; i++)
                {
                    Handles.color = (i % 3 == 0) ? anchorHandleColor : controlHandleColor;
                    Vector2 newPos = xyzToXZ(Handles.FreeMoveHandle(
                        xzToX0Z(spline[i]), Quaternion.identity, bezierHandleWidth, Vector3.zero, Handles.SphereHandleCap
                        ));
                    if (spline[i] != newPos)
                    {
                        Undo.RecordObject(this, "Move Path Point");
                        spline.MovePoint(i, newPos);
                        Debug.Log("Moving ctrl/anchor point " + i);
                        if (spline == roadSpline) TopologyHandler.OnSplineChanged();
                        if (OnUpdated != null) OnUpdated();
                    }
                }
            }
            #endregion
        }

        static Vector3 xzToX0Z(Vector2 xz)
        {
            return new Vector3(xz.x, 0, xz.y);
        }

        static Vector2 xyzToXZ(Vector3 xyz)
        {
            return new Vector2(xyz.x, xyz.z);
        }
    }

    EnergyCurve GeneratePolyline()
    {
        Debug.Log("Starting Road Generation...");

        ResetParameters();

        curve = new(Config);
        Debug.Log("Assign");
        return curve;
    }

    void ResetParameters()
    {
        numStuckIterations = 0;
        ls_currentStep = 0;
        subdivideCount = 0;

        steps = 0;
    }

    void ResetPolyline()
    {
        ResetParameters();
        curve.Reset(Config);
    }

    void PrintPolyline()
    {
        List<Vector2> polyline = curve.Polyline;
        string output = "";
        output += '[';
        for (int i = 0; i < polyline.Count; i++)
        {
            output += '[';
            output += polyline[i].x;
            output += ',';
            output += polyline[i].y;
            output += ']';
            if (i < polyline.Count - 1) output += ',';
        }
        output += ']';
        Debug.Log(output);
    }

    void StartEnvironment()
    {
        Debug.Log("Starting Environment...");
        scene = EditorSceneManager.NewScene(NewSceneSetup.DefaultGameObjects, NewSceneMode.Single);
        scene.name = "MapGen Environment";
    }

    void StopEnvironment()
    {
        EditorSceneManager.OpenScene("Assets/Scenes/Map4.unity");
    }

    void SaveObject()
    {
        Debug.Log("Saving Road Object...");
    }

    void RepulsionUpdate(int count)
    {
        for (int i = 0; i < count; i++)
            RepulsionUpdate();
    }

    void RepulsionUpdate()
    {
        steps++; 

        Tuple<bool, bool> returnedTuple = curve.ComputeLineSearchStep(repulsionType, useBarnesHut, useBackproj);
        bool goodStep = returnedTuple.Item1;
        bool shouldContinue = returnedTuple.Item2;
        if (!shouldContinue) runningLineSearch = false;

        // Only if we're actively running the simulation and not just doing a single step, because we want to be able to force more steps
        if (runningLineSearch)
        {
            ls_currentStep++;

            if (curve.normZero)
            {
                Debug.Log("Stopped because flow is near a local minimum.");
                runningLineSearch = false;
            }

            if (!goodStep)
            {
                numStuckIterations++;
                if (numStuckIterations >= 3 && curve.TargetLengthReached())
                {
                    Debug.Log("Stopped because flow hasn't made progress in a while");
                    runningLineSearch = false;
                }
            }
            else if (stepLimit > 0 && ls_currentStep >= stepLimit)
            {
                Debug.Log("Stopped because maximum number of steps was reached.");
                runningLineSearch = false;
            }
            else
            {
                numStuckIterations = 0;
            }
        }

        float avgLength = curve.TotalLength() / curve.NumEdges();
        if (avgLength > 2 * curve.initialAvgLength && subdivideCount < subdivideLimit)
        {
            //Debug.Log(avgLength);
            //Debug.Log(curve.initialAvgLength);
            subdivideCount++;
            SubdivideCurve();
        }

        Debug.Log("Step " + steps + " completed");
    }

    private void SubdivideCurve()
    {
        EnergyCurve subdivided = curve.Subdivide();
        //ReplaceNetwork(subdivided); //with network = subdivided;
    }

    void GeneateBezierSpline()
    {
        roadSpline = new(PolylineToBezier.Convert(curve.Polyline, epsilon, psi), true, true);
        //Debug.Log(bezierSpline);
    }

    void ScaleCurveToFitWidth()
    {
        float width = roadSpline.CalculateMinimumDistance();
        float newWidthMultiplier = width / RoadObjectWrapper.road.mesh.bounds.size.x;
        roadSpline.ScaleByT(WidthMultiplier / newWidthMultiplier);
    }

    void GenerateRoad()
    {
        RegenerateRoad();
    }

    void RegenerateRoad()
    {

        RoadGen.GenerateRoad_WholeParts(
            obj.transform,
            roadSpline,
            WidthMultiplier,
            HeightMultiplier,
            crossingShape,
            new()
            {
                {
                    TopologyType.Normal,
                    new(RoadObjectWrapper.road.mesh, RoadObjectWrapper.road.materials, roadPartLength) 
                },
                {
                    TopologyType.BridgeAscent,
                    new(RoadObjectWrapper.bridgeAscent.mesh, RoadObjectWrapper.bridgeAscent.materials, BridgeAscentLength())
                },
                {
                    TopologyType.Bridge,
                    new(RoadObjectWrapper.bridge.mesh, RoadObjectWrapper.bridge.materials, roadPartLength)
                },
                {
                    TopologyType.BridgeDescent,
                    new(RoadObjectWrapper.bridgeDescent.mesh, RoadObjectWrapper.bridgeDescent.materials, BridgeDescentLength())
                },
                {
                    TopologyType.Crossing,
                    new(RoadObjectWrapper.crossing.mesh, RoadObjectWrapper.crossing.materials, roadPartLength)
                },
                {
                    TopologyType.Ramp,
                    new(RoadObjectWrapper.ramp.mesh, RoadObjectWrapper.ramp.materials, RampLength())
                }
            }
        );

    }

    float RoadPartWidth()
    {
        return RoadObjectWrapper.road.mesh.bounds.size.x * WidthMultiplier;
    }

    float BridgeAscentLength()
    {
        return RoadObjectWrapper.bridgeAscent.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    float BridgeDescentLength()
    {
        return RoadObjectWrapper.bridgeAscent.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    float RampLength()
    {
        return RoadObjectWrapper.ramp.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    void LoadMap()
    {
        SerializedMapObject smo = SerializedMapObject.LoadFromFile(mapObjectName);
        roadSpline = smo.roadSpline;
        UpdateTopologyHandler(); // For new spline
        TopologyHandler.SetTopologies(smo.topologies, smo.fourPercCrossings_1D);
    }

    void SaveMap()
    {
        SerializedMapObject smo = new(roadSpline, TopologyHandler.topologies, TopologyHandler.fourPercCrossings);
        smo.SaveToFile(mapObjectName);
    }

    void ExportPrefab()
    {
        string baseObjectsPath = "Assets/Resources/Maps/Generated/Objects";
        AssetDatabase.CreateFolder(baseObjectsPath, mapObjectName);
        string prefabPath = $"{baseObjectsPath}/{mapObjectName}";
        AssetDatabase.CreateFolder(prefabPath, "Meshes");
        string meshesPath = $"{prefabPath}/Meshes";

        GameObject saveObject = Instantiate(RoadObjectWrapper.gameObject);
        saveObject.name = mapObjectName;
        saveObject.tag = "Untagged";
        DestroyImmediate(saveObject.transform.GetComponent<RoadObjectW>());

        var filters = saveObject.GetComponentsInChildren<MeshFilter>();
        for (int i = 0; i < filters.Length; i++)
        {
            MeshFilter filter = filters[i];
            Mesh mesh = filter.sharedMesh;
            AssetDatabase.CreateAsset(mesh, $"{meshesPath}/{i + 1}_{filter.gameObject.name}.asset"); // Creating the Asset suffices
        }

        PrefabUtility.SaveAsPrefabAsset(saveObject, $"{prefabPath}/{mapObjectName}.prefab", out bool success);
        if (success == true)
            Debug.Log("Prefab was saved successfully");
        else
            Debug.Log("Prefab failed to save" + success);

        DestroyImmediate(saveObject);
    }
}