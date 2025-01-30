using System;
using System.Collections.Generic;
using System.Security.Cryptography;
using UnityEditor;
using UnityEngine;
using UnityEngine.UIElements;

public class EditorMapGenTool : EditorWindow
{
    private MapGenTool mapGenTool;

    #region Drawing
    EnergyCurve Curve { get { return mapGenTool.GetEnergyCurve(); } }
    List<CurveVertex> PolyPoints { get { return Curve.verts; } }
    Spline Spline { get { return mapGenTool.roadSpline; } }
    #endregion

    /// <summary>
    /// When opening the window, add our Draw() function to the Unity drawing stack
    /// </summary>
    void OnEnable()
    {
        mapGenTool ??= new MapGenTool(); // only instantiate if not still available
        mapGenTool.InitAfterReload();

        SceneView.duringSceneGui += MyUpdateWrapper;
    }

    /// <summary>
    /// When closing the window, remove our Draw() function from the Unity drawing stack
    /// </summary>
    void OnDisable()
    {
        SceneView.duringSceneGui -= MyUpdateWrapper;
    }

    [MenuItem("Window/Map Gen Tool")]
    public static void ShowWindow()
    {
        GetWindow(typeof(EditorMapGenTool));
    }

    void OnGUI()
    {
        // You can only call GUI functions in OnGUI()
        mapGenTool.style_Section = new GUIStyle(EditorStyles.foldout)
        {
            fontStyle = FontStyle.Bold,
            alignment = TextAnchor.MiddleLeft,
            padding = new RectOffset(15, 0, 0, 0),
        };
        mapGenTool.style_Label = new GUIStyle(EditorStyles.label);

        mapGenTool.UpdateUI();
    }

    void MyUpdateWrapper(SceneView sceneView)
    {
        if (mapGenTool.genModeConfig is not GenModeConfig_Bezier && mapGenTool.genModeConfig is not GenModeConfig_Circular)
        {
            Debug.Log("Tool not yet initialized, retrying...");
            OnGUI(); // OnGUI necessary to initialize some things. Only called when window is visible
            return;
        }

        try
        {
            mapGenTool.MyUpdate();
            Draw();
        }
        catch (NullReferenceException)
        {
            Debug.Log("No curve found, generating...");
            mapGenTool.GeneratePolyline();
            Debug.Log("Curve generated");
        }

        mapGenTool.GetEnergyCurve().SerializeVertsEdgesPositions();
        mapGenTool.GetEnergyCurve().S_SerializeObstaclePositions();
    }

    #region Drawing
    void Draw()
    {
        #region Colors

        Color controlHandleColor = new(0.2f, 0.5f, 0.2f);
        Color anchorHandleColor = new(255f / 255f, 225f / 255f, 67f / 255f);
        Color lineColor = new(0, 0, 0);
        Color splineColor = new(0.2f, 0.2f, 0.5f);
        Color polyColor = new(0.4f, 0.6f, 1f);

        float bezierWidth = 8f; // 20f
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
            if (mapGenTool.showPolyLine)
            {
                List<Vector3> polyPoints3D = new List<Vector3>();
                for (int i = 0; i < PolyPoints.Count; i++)
                    polyPoints3D.Add(xzToX0Z(PolyPoints[i].Position()));
                if (mapGenTool.curveClosed) polyPoints3D.Add(xzToX0Z(PolyPoints[0].Position())); // To draw a closed curve, we need to add the first point again
                Handles.DrawAAPolyLine(EditorGUIUtility.whiteTexture, polylineWidth, polyPoints3D.ToArray());
            }
            if (mapGenTool.showPolyPoints)
            {
                for (int i = 0; i < PolyPoints.Count; i++)
                {
                    Handles.DrawSolidDisc(xzToX0Z(PolyPoints[i].Position()), Vector3.up, polylineDiscWidth);
                }
            }
        }
        #endregion

        #region Obstacles
        if (mapGenTool.showObstacles && Curve != null && Curve.obstacles != null)
        {
            Handles.color = Color.red;
            foreach (Obstacle obs in Curve.obstacles)
            {
                if (!obs.IsEnabled)
                    continue;

                List<Vector3> obstaclePoints3D = new List<Vector3>();
                foreach (CurveVertex vert in obs.verts)
                    obstaclePoints3D.Add(xzToX0Z(vert.Position()));
                obstaclePoints3D.Add(xzToX0Z(obs.verts[0].Position())); // To draw a closed curve, we need to add the first point again
                if (mapGenTool.showPolyLine)
                {
                    Handles.DrawAAPolyLine(EditorGUIUtility.whiteTexture, obstacleWidth, obstaclePoints3D.ToArray());
                }
                if (mapGenTool.showPolyPoints)
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
        if (Spline != null && mapGenTool.showBezier)
        {
            DrawSpline(Spline);
        }
        #endregion

        void DrawSpline(Spline spline)
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

                if (mapGenTool.showBezierHandles)
                {
                    Handles.DrawLine(points3D[1], points3D[0], bezierLineWidth);
                    Handles.DrawLine(points3D[2], points3D[3], bezierLineWidth);
                }
                Handles.DrawBezier(points3D[0], points3D[3], points3D[1], points3D[2], splineColor, Texture2D.whiteTexture, bezierWidth);
            }
            #endregion

            #region Handles
            if (mapGenTool.showBezierHandles)
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
                        if (spline == Spline) TopologyHandler.OnSplineChanged();
                        if (mapGenTool.autoUpdateRoad) mapGenTool.RegenerateRoad();
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
    #endregion

    [Serializable] // Serializable here to save and keep everything (what we can) when editor reloads / scripts recompile
    private class MapGenTool : MapGenToolBase
    {
        Vector2 scrollPos;

        public GUIStyle style_Section;
        public GUIStyle style_Label;

        Dictionary<int, bool> foldouts;

        GameObject obj; // As we can't serialize, we store everything in an object in the scene

        RoadObjectW _roadObjectWrapper;
        RoadObjectW RoadObjectWrapper // MonoBehavior for forced Serialization
        {
            get
            {
                if (_roadObjectWrapper == null)
                    _roadObjectWrapper = obj.GetComponent<RoadObjectW>();
                return _roadObjectWrapper;
            }
        }
        protected override MeshTopologyEditorConfig GetObject_road() => RoadObjectWrapper.road;
        protected override MeshTopologyEditorConfig GetObject_bridgeAscent() => RoadObjectWrapper.bridgeAscent;
        protected override MeshTopologyEditorConfig GetObject_bridge() => RoadObjectWrapper.bridge;
        protected override MeshTopologyEditorConfig GetObject_bridgeDescent() => RoadObjectWrapper.bridgeDescent;
        protected override MeshTopologyEditorConfig GetObject_crossing() => RoadObjectWrapper.crossing;
        protected override MeshTopologyEditorConfig GetObject_ramp() => RoadObjectWrapper.ramp;

        public override EnergyCurve GetEnergyCurve() {
            if (obj == null)
                InitObj();
            return obj.GetComponent<EnergyCurveW>().curve;
        }
        protected override void SetEnergyCurve(EnergyCurve curve) => obj.GetComponent<EnergyCurveW>().curve = curve;

        public MapGenTool()
        {
            foldouts = new();
        }

        public void InitAfterReload()
        {
            try
            {
                GetEnergyCurve().TryInitVertsEdgesFromSerialized();
                GetEnergyCurve().S_TryInitObstacles();
            }
            catch (Exception)
            {
                // Wait for MyUpdate to generate polyline
            }
        }

        void InitObj()
        {
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

        public override void UpdateUI()
        {
            GUILayout.BeginVertical();
            GUILayout.Space(10);
            scrollPos = EditorGUILayout.BeginScrollView(scrollPos);

            base.UpdateUI();

            EditorGUILayout.EndScrollView();
            GUILayout.EndVertical();
        }

        public override void RegenerateRoad()
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

        protected override void LoadMap()
        {
            SerializedMapObject smo = SerializedMapObject.LoadFromFile(mapObjectName);
            roadSpline = smo.roadSpline;
            UpdateTopologyHandler(); // For new spline
            TopologyHandler.SetTopologies(smo.topologies, smo.fourPercCrossings_1D);
        }

        protected override void SaveMap()
        {
            SerializedMapObject smo = new(roadSpline, TopologyHandler.topologies, TopologyHandler.fourPercCrossings);
            smo.SaveToFile(mapObjectName);
        }

        protected override void ExportPrefab()
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

        #region GUI methods
        protected override void CreateSection(string name, Action content)
        {
            CreateFoldout(name, content, 20, 6);
        }

        protected override void CreateSubSection(string name, Action content)
        {
            CreateFoldout(name, content, 16, 4);
        }

        protected override void CreateSubSubSection(string name, Action content)
        {
            CreateFoldout(name, content, 14, 3);
        }

        protected override void CreateLabel(string name, int fontSize = 12)
        {
            GUIStyle style = new(style_Label) { fontSize = fontSize };
            GUILayout.Label(name, style);
        }

        protected override void CreateTextField(string name, ref string input)
        {
            input = EditorGUILayout.TextField(name, input);
        }

        protected override void CreateIntField(string name, ref int input)
        {
            input = EditorGUILayout.IntField(name, input);
        }

        protected override void CreateFloatField(string name, ref float input)
        {
            input = EditorGUILayout.FloatField(name, input);
        }

        protected override void CreateVector2Field(string name, ref Vector2 input)
        {
            input = EditorGUILayout.Vector2Field(name, input);
        }

        protected override void CreateButton(string name, Action callback)
        {
            if (GUILayout.Button(name)) callback();
        }

        protected override void CreateCheckbox(string name, ref bool trigger)
        {
            trigger = GUILayout.Toggle(trigger, name);
        }

        protected override void CreateCheckbox_Dict<T>(string name, Dictionary<T, bool> dict, T key)
        {
            dict[key] = GUILayout.Toggle(dict[key], name);
        }

        protected override void CreateSlider(ref int value, string text, int start, int end)
        {
            Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
            position.height *= 0.5f;
            value = (int)EditorGUI.Slider(position, text, value, start, end);

            position.y += position.height;
            position.x += EditorGUIUtility.labelWidth;
            position.width -= EditorGUIUtility.labelWidth + 54;

            GUIStyle style = GUI.skin.label;
            style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
            style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
        }

        protected override void CreateSliderFloat(ref float value, string text, float start, float end, float snap = 0.5f)
        {
            Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
            position.height *= 0.5f;
            value = EditorGUI.Slider(position, text, value, start, end);

            position.y += position.height;
            position.x += EditorGUIUtility.labelWidth;
            position.width -= EditorGUIUtility.labelWidth + 54;

            GUIStyle style = GUI.skin.label;
            style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
            style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
        }

        protected override void CreateCurveField(string name, ref AnimationCurve curve)
        {
            curve = EditorGUILayout.CurveField(name, curve);
        }

        protected override void CreateFoldout(string name, Action content, int fontSize, int bottomSpacing = 0)
        {
            int foldoutHash = content.GetHashCode();

            GUIStyle style = new(style_Section) { fontSize = fontSize };

            bool initialization = false;
            if (!foldouts.ContainsKey(foldoutHash))
            {
                initialization = true;
                foldouts.Add(foldoutHash, true); // Initially everything must be folded out, as otherwise some initializations fail
            }

            foldouts[foldoutHash] = EditorGUILayout.Foldout(
                foldout: foldouts.ContainsKey(foldoutHash) && foldouts[foldoutHash],
                content: name,
                toggleOnLabelClick: true,
                style: style
            );

            if (foldouts[foldoutHash])
            {
                GUILayout.BeginHorizontal();
                GUILayout.Space(10);

                GUILayout.BeginVertical();
                GUILayout.Space(bottomSpacing);

                content();

                GUILayout.Space(bottomSpacing); // 2. Additional spacing under last element
                GUILayout.EndVertical();
                GUILayout.EndHorizontal();
            }
            GUILayout.Space(bottomSpacing); // 1. Default bottom spacing when collapsed

            if (initialization)
            {
                foldouts[foldoutHash] = false;
            }
        }

        protected override void CreateEnumSelection<T>(string name, T currentSelection, Action<T> onValueChanged)
        {
            currentSelection = (T)EditorGUILayout.EnumPopup("Gen Mode", currentSelection);
            onValueChanged(currentSelection);
        }
        #endregion
    }
}
