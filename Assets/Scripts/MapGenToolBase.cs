using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using UnityEditor;
using UnityEngine;

public abstract class MapGenToolBase
{
    #region Automation Controls
    Auto_Stage auto_stage = Auto_Stage.PolylineGen;
    int auto_currCurveCount = 0;
    int auto_intersectionCount = 0;

    #region Exposed
    int auto_maxCurveCount = 5;
    bool auto_generating = false;
    string auto_filepath = "C:\\temp\\";
    #endregion

    #endregion

    #region Curve

    #region Polyline

    #region Exposed
    public GenMode genMode = GenMode.Circular;
    public RepulsionType repulsionType = RepulsionType.Sobolev;

    public bool curveClosed = true;
    bool deacObsAfterScaling = true;
    bool rotateAfterScaling = false;
    bool noRepulsionAfterScaling = false;
    float lsStepThreshold = 1e-15f;
    float energyThreshold = 0f;
    public bool showBezier = true;
    public bool showBezierHandles = true;
    public bool showPolyLine = true;
    public bool showPolyPoints = true;
    public bool showObstacles = true;
    bool useBarnesHut = true;
    bool useBackproj = true;
    bool runningLineSearch = false;

    [SerializeField] // Serialization is important for editor tool
    public GenModeConfig genModeConfig = new GenModeConfig_Circular();

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

    public abstract EnergyCurve GetEnergyCurve();
    protected abstract void SetEnergyCurve(EnergyCurve curve);
    #endregion

    #region BezierSpline
    // Spline
    public RoadSpline roadSpline;
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
    protected abstract MeshTopologyEditorConfig GetObject_road();
    protected abstract MeshTopologyEditorConfig GetObject_bridgeAscent();
    protected abstract MeshTopologyEditorConfig GetObject_bridge();
    protected abstract MeshTopologyEditorConfig GetObject_bridgeDescent();
    protected abstract MeshTopologyEditorConfig GetObject_crossing();
    protected abstract MeshTopologyEditorConfig GetObject_ramp();

    protected float roadPartLength = 1f;
    float bridgePartLength = 1f;
    float _widthMultiplier = 1f;
    protected float WidthMultiplier { get { return _widthMultiplier * 1000; } }
    float _heightMultiplier = 1f;
    protected float HeightMultiplier { get { return _heightMultiplier * 1000; } }
    protected float crossingShape = 0.552f;

    public bool autoUpdateRoad = true;

    protected string mapObjectName = "MyMap";
    #endregion

    public virtual void UpdateUI()
    {
        #region Automation
        CreateSection("Automation", () =>
        {
            CreateIntField("Number", ref auto_maxCurveCount);
            CreateLabel($"Already Generated: {auto_currCurveCount}");
            CreateButton("Generate", Auto_Start);
            CreateButton("Stop", Auto_Stop);
            CreateTextField("Export Directory", ref auto_filepath);
        });
        #endregion

        GUILayout.Space(10);

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
            CreateButton("Add Intersection", () => AddIntersection());

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
                    CreateEnumSelection("Gen Mode", genMode, (GenMode value) => genMode = value);
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

                        //EditorGUIUtility.labelWidth = 50;
                        //EditorGUIUtility.fieldWidth = 20;
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
                        //EditorGUIUtility.labelWidth = 0;
                        //EditorGUIUtility.fieldWidth = 0;
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
                    //EditorGUIUtility.labelWidth = 50;
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
                    //EditorGUIUtility.labelWidth = 0;

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
                    CreateEnumSelection("Repulsion Method", repulsionType, (RepulsionType value) => repulsionType = value);
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
#if UNITY_EDITOR
                // We use Asset selections via EditorGUILayout.ObjectField(),
                // which is only possible in the Unity Editor.
                // Note that this doesn't really work for the Play Mode Editor,
                // just for the Unity Editor Editor. For the Play Mode Editor,
                // you should assign the assets in the MonoBehavior directly.

                // Note that the "name" arguments in the DoFoldout() calls also ensures
                // the hashing in CreateSubSubSection to work correctly

                CreateSubSubSection("Road Config", () =>
                {
                    DoFoldout("Road", GetObject_road());
                });

                CreateSubSubSection("BridgeAscent Config", () =>
                {
                    DoFoldout("BridgeAscent", GetObject_bridgeAscent());
                });

                CreateSubSubSection("Bridge Config", () =>
                {
                    DoFoldout("Bridge", GetObject_bridge());
                });

                CreateSubSubSection("BridgeDescent Config", () =>
                {
                    DoFoldout("BridgeDescent", GetObject_bridgeDescent());
                });

                CreateSubSubSection("Crossing Config", () =>
                {
                    DoFoldout("Crossing", GetObject_crossing());
                });

                CreateSubSubSection("Ramp Config", () =>
                {
                    DoFoldout("Ramp", GetObject_ramp());
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

#else
                CreateLabel("Assigning objects here is not supported yet!");
                CreateLabel("Please assign road objects in the Unity Editor");
#endif

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

        UpdateTopologyHandler();
    }

    public void MyUpdate()
    {
        if (auto_generating)
        {
            switch (auto_stage)
            {
                case Auto_Stage.PolylineGen:
                    auto_intersectionCount = 0;
                    GeneratePolyline();

                    runningLineSearch = true;
                    auto_stage = Auto_Stage.LineSearch;
                    break;
                case Auto_Stage.LineSearch:
                    if (!runningLineSearch) // generation done
                    {
                        auto_stage = Auto_Stage.SplineGen;
                    }
                    break;
                case Auto_Stage.SplineGen:
                    GeneateBezierSpline();
                    auto_stage = Auto_Stage.IntersectionGen;
                    break;
                case Auto_Stage.IntersectionGen:
                    bool intersectionAdded = AddIntersection();
                    if (intersectionAdded)
                        auto_intersectionCount++;
                    else
                        auto_stage = Auto_Stage.FittingWidth;
                    break;
                case Auto_Stage.FittingWidth:
                    ScaleCurveToFitWidth();
                    auto_stage = Auto_Stage.Finished;
                    break;
                case Auto_Stage.Finished:
                    ExportTrack();
                    auto_currCurveCount++;
                    if (auto_currCurveCount < auto_maxCurveCount)
                        auto_stage = Auto_Stage.PolylineGen;
                    else
                        auto_generating = false;
                    break;
                default:
                    Debug.LogWarning("Something went wrong!");
                    break;
            }
        }

        if (runningLineSearch)
        {
            RepulsionUpdate();
        }

        GetEnergyCurve().SerializeVertsEdgesPositions();
        GetEnergyCurve().S_SerializeObstaclePositions();
    }

    protected void UpdateTopologyHandler()
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

    private void Auto_Start()
    {
        auto_currCurveCount = 0;
        if (auto_currCurveCount < auto_maxCurveCount)
        {
            auto_generating = true;
            auto_stage = Auto_Stage.PolylineGen;
        }
    }

    private void Auto_Stop()
    {
        auto_generating = false;
        runningLineSearch = false;
    }

    private void ExportTrack()
    {
        string timestamp = DateTime.Now.ToString("yyyyMMddHHmmss");

        string points = SplineToPointsString();
        string pointsFilename = $"{timestamp}_{lengthScale}_{deacObsAfterScaling}_points.csv";
        string csv_points = GenerateCsvContent(pointsFilename, points);

        Debug.Log(csv_points.ToString());
        try
        {
            File.WriteAllText(auto_filepath + "\\" + pointsFilename, csv_points.ToString());
        }
        catch (Exception ex)
        {
            Debug.LogError("Exporting points failed!");
            Debug.LogException(ex);
        }


        string controls = SplineToControlsString();
        string controlsFilename = $"{timestamp}_{lengthScale}_{deacObsAfterScaling}_spline.csv";
        string csv_controls = GenerateCsvContent(controlsFilename, controls);

        Debug.Log(csv_controls);
        try
        {
            File.WriteAllText(auto_filepath + "\\" + controlsFilename, csv_controls.ToString());
        }
        catch (Exception ex)
        {
            Debug.LogError("Exporting controls failed!");
            Debug.LogException(ex);
        }
    }

    private string GenerateCsvContent(string filename, string data)
    {
        return string.Join(";", new[]
        {
            filename,
            DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"),
            steps.ToString(),
            lengthScale.ToString(), // "complexity"
            roadSpline.TotalLength.ToString(), // actual length
            deacObsAfterScaling.ToString(), // isometry (at least we assume so, but could be not related)
            auto_intersectionCount.ToString(), // number of deliberate intersections
            roadSpline.CalculateIntersections().Count.ToString(), // number if actual intersections
            data
        });
    }

    private string SplineToPointsString()
    {
        StringBuilder sb = new();
        bool first = true;
        foreach (Vector2 v in roadSpline.CalculateEvenlySpacedPoints(1.0f))
        {
            if (!first)
                sb.Append(", ");
            first = false;
            //Don't us Vector2 toString because it rounds to two decimals
            sb.Append("(" + v.x.ToString(CultureInfo.InvariantCulture) + ", " + v.y.ToString(CultureInfo.InvariantCulture) + ")");
        }
        return sb.ToString();
    }

    private string SplineToControlsString()
    {
        StringBuilder sb = new();
        bool first = true;
        foreach (Vector2[] segControls in roadSpline.pointsInSegment)
        {
            foreach (Vector2 c in segControls) // four controls per segment -> anchor points are always included in two segements
            {
                if (!first)
                    sb.Append(", ");
                first = false;
                sb.Append("(" + c.x.ToString(CultureInfo.InvariantCulture) + ", " + c.y.ToString(CultureInfo.InvariantCulture) + ")");
            }
        }
        return sb.ToString();
    }

    public EnergyCurve GeneratePolyline()
    {
        ResetParameters();

        SetEnergyCurve(new(Config));
        return GetEnergyCurve();
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
        GetEnergyCurve().Reset(Config);
    }

    void PrintPolyline()
    {
        List<Vector2> polyline = GetEnergyCurve().Polyline;
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

    void RepulsionUpdate(int count)
    {
        for (int i = 0; i < count; i++)
            RepulsionUpdate();
    }

    void RepulsionUpdate()
    {
        steps++;

        (bool goodStep, bool shouldContinue) = GetEnergyCurve().ComputeLineSearchStep(repulsionType, useBarnesHut, useBackproj);
        if (!shouldContinue) runningLineSearch = false;

        // Only if we're actively running the simulation and not just doing a single step, because we want to be able to force more steps
        if (runningLineSearch)
        {
            ls_currentStep++;

            if (GetEnergyCurve().normZero)
            {
                Debug.Log("Stopped because flow is near a local minimum.");
                runningLineSearch = false;
            }

            if (!goodStep)
            {
                numStuckIterations++;
                if (numStuckIterations >= 3 && GetEnergyCurve().TargetLengthReached())
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

        float avgLength = GetEnergyCurve().TotalLength() / GetEnergyCurve().NumEdges();
        if (avgLength > 2 * GetEnergyCurve().initialAvgLength && subdivideCount < subdivideLimit)
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
        EnergyCurve subdivided = GetEnergyCurve().Subdivide();
        //ReplaceNetwork(subdivided); //with network = subdivided;
    }

    void GeneateBezierSpline()
    {
        roadSpline = new(PolylineToBezier.Convert(GetEnergyCurve().Polyline, epsilon, psi), true, true);
        //Debug.Log(bezierSpline);
    }

    bool AddIntersection()
    {
        return roadSpline.AddIntersections(intersectionPreferredDistance);
    }

    void ScaleCurveToFitWidth()
    {
        float width = roadSpline.CalculateMinimumDistance();
        float newWidthMultiplier = width / GetObject_road().mesh.bounds.size.x;
        roadSpline.ScaleByT(WidthMultiplier / newWidthMultiplier);
    }

    void GenerateRoad()
    {
        RegenerateRoad();
    }

    float RoadPartWidth()
    {
        return GetObject_road().mesh.bounds.size.x * WidthMultiplier;
    }

    protected float BridgeAscentLength()
    {
        return GetObject_bridgeAscent().mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    protected float BridgeDescentLength()
    {
        return GetObject_bridgeAscent().mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    protected float RampLength()
    {
        return GetObject_ramp().mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    #region abstract methods
    public abstract void RegenerateRoad();
    protected abstract void LoadMap();
    protected abstract void SaveMap();
    protected abstract void ExportPrefab();

    #region GUI methods
    protected abstract void CreateSection(string name, Action content);
    protected abstract void CreateSubSection(string name, Action content);
    protected abstract void CreateSubSubSection(string name, Action content);
    protected abstract void CreateLabel(string name, int fontSize = 12);
    protected abstract void CreateTextField(string name, ref string input);
    protected abstract void CreateIntField(string name, ref int input);
    protected abstract void CreateFloatField(string name, ref float input);
    protected abstract void CreateVector2Field(string name, ref Vector2 input);
    protected abstract void CreateButton(string name, Action callback);
    protected abstract void CreateCheckbox(string name, ref bool trigger);
    protected abstract void CreateCheckbox_Dict<T>(string name, Dictionary<T, bool> dict, T key);
    protected abstract void CreateSlider(ref int value, string text, int start, int end);
    protected abstract void CreateSliderFloat(ref float value, string text, float start, float end, float snap = 0.5f);
    protected abstract void CreateCurveField(string name, ref AnimationCurve curve);
    protected abstract void CreateFoldout(string name, Action content, int fontSize, int bottomSpacing = 0);
    protected abstract void CreateEnumSelection<T>(string name, T currentSelection, Action<T> onValueChanged) where T : Enum;
    #endregion
    #endregion
}

enum Auto_Stage
{
    PolylineGen,
    LineSearch,
    SplineGen,
    IntersectionGen,
    FittingWidth,
    Finished
}