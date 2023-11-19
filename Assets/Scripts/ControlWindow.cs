using UnityEngine;
using System.Collections.Generic;
using System;
using System.Linq;

public class ControlWindow : ToolWindow
{
    private Vector2 scrollPosition;

    #region Curve

    #region Polyline

    #region Exposed

    public GenMode genMode = GenMode.Circular;
    public RepulsionType repulsionType = RepulsionType.Normal;

    public bool curveClosed = true;
    bool deacObsAfterScaling = true;
    bool rotateAfterScaling = false;
    bool noRepulsionAfterScaling = false;
    bool showBezier = true;
    bool showBezierHandles = true;
    public bool showPolyLine = true;
    public bool showPolyPoints = true;
    bool showObstacles = true;
    bool useBarnesHut = true;
    bool useBackproj = true;
    bool runningLineSearch = false;

    [SerializeField]
    GenModeConfig genModeConfig = new GenModeConfig_Circular();

    public EnergyCurve_EditorConfig Config
    {
        get
        {
            return new EnergyCurve_EditorConfig(curveClosed, lengthScale, deacObsAfterScaling, rotateAfterScaling, noRepulsionAfterScaling, genModeConfig, ObstacleList, numObstacles, PotentialList, ConstraintList);
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

    EnergyCurve _curve;
    [HideInInspector]
    public EnergyCurve Curve
    {
        get
        {
            return _curve ??= GeneratePolyline();
        }
        set
        {
            _curve = value;
        }
    }
    Spline InitialSpline { get { return Curve.spline; } }

    #endregion

    #region BezierSpline

    // Spline
    RoadSpline roadSpline;
    float epsilon = 0.2f;
    float psi = 2f;

    // Topology
    float _fixedPartLengthMultiplier = 1f;
    float FixedPartLengthMultiplier { get { return _fixedPartLengthMultiplier * 1000; } }
    float crossingExtraSize = 0.2f;

    #endregion

    #endregion

    #region RoadGen

    bool overrideWidth = false;
    float roadPartLength = 1f;
    float bridgePartLength = 1f;
    float _widthMultiplier = 1f;
    float WidthMultiplier { get { return _widthMultiplier * 1000; } }
    float _heightMultiplier = 8f;
    float HeightMultiplier { get { return _heightMultiplier * 1000; } }
    float crossingShape = 0.5f;

    bool autoUpdateRoad = true;

    public MeshTopologyEditorConfig object_road;
    public MeshTopologyEditorConfig object_bridgeAscent;
    public MeshTopologyEditorConfig object_bridge;
    public MeshTopologyEditorConfig object_bridgeDescent;
    public MeshTopologyEditorConfig object_crossing;
    public MeshTopologyEditorConfig object_ramp;

    #endregion

    /*
    [MenuItem("Window/Map Gen -- Default Values")]
    public static void ShowWindow()
    {
        GetWindow(typeof(DefaultValues));
    }
    */

    protected override void OnGUI()
    {
        base.OnGUI();

        float maxWidth = 300;
        float contentHeight = 2000;

        // Background
        GUI.Box(new Rect(0, 0, maxWidth, Screen.height), GUIContent.none);

        GUILayout.BeginVertical();

        // Begin the ScrollView
        scrollPosition = GUI.BeginScrollView(
            new Rect(0, 0, maxWidth, Screen.height),
            scrollPosition,
            new Rect(0, 0, maxWidth - 20, contentHeight) // Subtracting 20 for the scrollbar width
        );

        #region Options

        CreateSection("Options", () =>
        {

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
                    //Debug.LogWarning("GenMode Dropdown tbd");
                    //genMode = (GenMode)EditorGUILayout.EnumPopup("Gen Mode", genMode);
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
                            CreateVector3Field("Center", ref obs.center);

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

                    #endregion

                    #region Repulsion
                    CreateLabel("");
                    CreateLabel("Repulsion");
                    CreateCheckbox("Use Barnes Hut", ref useBarnesHut);
                    CreateCheckbox("Use Backprojection", ref useBackproj);
                    //repulsionType = (RepulsionType)EditorGUILayout.EnumPopup("Repulsion Method", repulsionType);
                    //Debug.LogWarning("RepulsionType Dropdown tbd");

                    #endregion

                });

                #endregion

                #region Road Spline
                CreateSubSubSection("Bezier Spline", () =>
                {
                    CreateFloatField("Epsilon", ref epsilon);
                    CreateFloatField("Psi", ref psi);
                    CreateFloatField("Fixed Part Length Multiplier", ref _fixedPartLengthMultiplier);
                    CreateFloatField("Crossing Extra Size", ref crossingExtraSize);
                });

                #endregion
            });
            #region Road

            CreateSubSection("Road Object", () =>
            {
                DoFoldout("Road", object_road);
                DoFoldout("BridgeAscent", object_bridgeAscent);
                DoFoldout("Bridge", object_bridge);
                DoFoldout("BridgeDescent", object_bridgeDescent);
                DoFoldout("Crossing", object_crossing);
                DoFoldout("Ramp", object_ramp);

                void DoFoldout(string name, MeshTopologyEditorConfig metc)
                {
                    metc ??= new(null, 0, new(), false);
                    _CreateFoldout(name + " Config", ref metc.showFoldout);
                    if (metc.showFoldout)
                    {
                        //Debug.LogWarning("Mesh Dropdown tbd");
                        //metc.mesh = (Mesh)EditorGUILayout.ObjectField(name + " Mesh", metc.mesh, typeof(Mesh), allowSceneObjects: false);
                        CreateSlider(ref metc.numMaterials, name + " Materials", 0, 5);
                        CreateMaterialSelection(metc.materials, metc.numMaterials);
                    }
                }

                void CreateMaterialSelection(List<Material> matList, int numMats)
                {
                    for (int i = 0; i < numMats; i++)
                    {
                        if (i >= matList.Count) matList.Add(new(Shader.Find("Specular")));
                        //Debug.LogWarning("Material Dropdown tbd");
                        //matList[i] = (Material)EditorGUILayout.ObjectField("Material " + (i + 1), matList[i], typeof(Material), allowSceneObjects: false);
                    }

                    if (matList != null && matList.Count > numMats)
                    {
                        matList.RemoveRange(Mathf.Max(numMats - 1, 0), matList.Count - numMats);
                    }
                }

                CreateCheckbox("Override Width", ref overrideWidth);
                if (overrideWidth) CreateFloatField("Width Multiplier", ref _widthMultiplier);
                CreateFloatField("Height Multiplayer", ref _heightMultiplier);
                CreateFloatField("Road Part Length", ref roadPartLength);
                CreateFloatField("Bridge Part Length", ref bridgePartLength);
                CreateSliderFloat(ref crossingShape, "Crossing Shape", 0f, 1f);
                CreateCheckbox("Auto Update Road", ref autoUpdateRoad);

            });

            #endregion



        });

        #endregion

        #region Running

        CreateSection("Running", () =>
        {
            CreateButton("Create Polyline", () => GeneratePolyline());
            CreateButton("Reset Polyline", ResetPolyline);
            CreateButton("Print Polyline as 2D Arary", PrintPolyline);

            CreateButton("Single Step LS", () => RepulsionUpdate());
            CreateButton("Ten Steps LS", () => RepulsionUpdate(10));
            CreateCheckbox("Run LS", ref runningLineSearch);
        });

        #endregion

        GUI.EndScrollView();

        GUILayout.EndVertical();
    }


    void Update()
    {
        if (runningLineSearch)
        {
            RepulsionUpdate();
        }

        Curve.SerializeVertsEdgesPositions();
        Curve.S_SerializeObstaclePositions();
    }

    EnergyCurve GeneratePolyline()
    {
        Debug.Log("Generating Polyline");
        ResetParameters();

        Curve = new(Config);

        return Curve;
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
        Curve.Reset(Config);
    }

    void PrintPolyline()
    {
        List<Vector3> polyline = Curve.Polyline;
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

        Tuple<bool, bool> returnedTuple = Curve.ComputeLineSearchStep(repulsionType, useBarnesHut, useBackproj);
        bool goodStep = returnedTuple.Item1;
        bool shouldContinue = returnedTuple.Item2;
        if (!shouldContinue) runningLineSearch = false;

        // Only if we're actively running the simulation and not just doing a single step, because we want to be able to force more steps
        if (runningLineSearch)
        {
            ls_currentStep++;

            if (Curve.normZero)
            {
                Debug.Log("Stopped because flow is near a local minimum.");
                runningLineSearch = false;
            }

            if (!goodStep)
            {
                numStuckIterations++;
                if (numStuckIterations >= 3 && Curve.TargetLengthReached())
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

        float avgLength = Curve.TotalLength() / Curve.NumEdges();
        if (avgLength > 2 * Curve.initialAvgLength && subdivideCount < subdivideLimit)
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
        EnergyCurve subdivided = Curve.Subdivide();
        //ReplaceNetwork(subdivided); //with network = subdivided;
    }

    void GeneateBezierSpline()
    {
        roadSpline = new(PolylineToBezier.Convert(Curve.Polyline, epsilon, psi), true, true);
        //Debug.Log(bezierSpline);
    }

    void CalculateWidthMultiplier()
    {
        float width = roadSpline.CalculateMinimumDistance();
        float new_widthMultiplier = width / object_road.mesh.bounds.size.x;
        float _wToW = WidthMultiplier / _widthMultiplier;
        _widthMultiplier = new_widthMultiplier / _wToW;
        Debug.Log("New WidthMultiplier: " + WidthMultiplier);
    }

    void RegenerateRoad()
    {
        /*
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
        */

    }

    float RoadPartWidth()
    {
        return object_road.mesh.bounds.size.x * WidthMultiplier;
    }

    float BridgeAscentLength()
    {
        return object_bridgeAscent.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    float BridgeDescentLength()
    {
        return object_bridgeAscent.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    float RampLength()
    {
        return object_ramp.mesh.bounds.size.z * FixedPartLengthMultiplier;
    }

    /*

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

    */
}

[Serializable] public class Dict_ConstraintType_Bool : SerializableDictionary<ConstraintType, bool> { }
[Serializable] public class Dict_PotentialType_Bool : SerializableDictionary<PotentialType, bool> { }
[Serializable] public class Dict_PotentialType_PotentialConfig : SerializableDictionary<PotentialType, PotentialConfig> { }

[Serializable]
public class SerializableDictionary<TKey, TValue> : Dictionary<TKey, TValue>, ISerializationCallbackReceiver
{
    [SerializeField]
    private List<TKey> keys = new List<TKey>();

    [SerializeField]
    private List<TValue> values = new List<TValue>();

    // save the dictionary to lists
    public void OnBeforeSerialize()
    {
        keys.Clear();
        values.Clear();
        foreach (KeyValuePair<TKey, TValue> pair in this)
        {
            keys.Add(pair.Key);
            values.Add(pair.Value);
        }
    }

    // load dictionary from lists
    public void OnAfterDeserialize()
    {
        this.Clear();

        if (keys.Count != values.Count)
            throw new System.Exception(string.Format("there are {0} keys and {1} values after deserialization. Make sure that both key and value types are serializable."));

        for (int i = 0; i < keys.Count; i++)
            this.Add(keys[i], values[i]);
    }
}

[Serializable]
public class MeshTopologyEditorConfig
{
    public Mesh mesh;
    public int numMaterials = 3;
    public List<Material> materials = new();
    public bool showFoldout;

    public MeshTopologyEditorConfig(Mesh mesh, int numMaterials, List<Material> materials, bool showFoldout)
    {
        this.mesh = mesh;
        this.numMaterials = numMaterials;
        this.materials = materials;
        this.showFoldout = showFoldout;
    }
}
