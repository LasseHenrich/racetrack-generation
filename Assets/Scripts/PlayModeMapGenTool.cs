using System;
using System.Collections.Generic;
using UnityEngine;

public class PlayModeMapGenTool : MonoBehaviour
{
    private MapGenTool mapGenTool;

    [SerializeField] MeshTopologyEditorConfig object_road;
    [SerializeField] MeshTopologyEditorConfig object_bridgeAscent;
    [SerializeField] MeshTopologyEditorConfig object_bridge;
    [SerializeField] MeshTopologyEditorConfig object_bridgeDescent;
    [SerializeField] MeshTopologyEditorConfig object_crossing;
    [SerializeField] MeshTopologyEditorConfig object_ramp;

    #region Drawing
    [SerializeField] Material mat;

    #region Colors

    [SerializeField] Color controlHandleColor = new(0.2f, 0.5f, 0.2f);
    [SerializeField] Color anchorHandleColor = new(1.0f, 1.0f, 67f / 255f);
    [SerializeField] Color lineColor = new(0, 0, 0);
    [SerializeField] Color splineColor = new(0.2f, 0.2f, 0.5f);
    [SerializeField] Color polyColor = new(0.4f, 0.6f, 1f);
    [SerializeField] Color obstacleColor = new(1.0f, 0f, 0f);

    #endregion

    #region Sizes

    [SerializeField] float bezierWidth = 8f; // 20f
    [SerializeField] float bezierLineWidth = 3f; // 5f
    [SerializeField] float bezierHandleWidth = 1.4f;//1.0f;//0.6f;
    [SerializeField] float polylineWidth = 4f;
    [SerializeField] float polylineDiscWidth = 0.15f;
    [SerializeField] float obstacleWidth = 2f;
    [SerializeField] float obstacleDiscWidth = 0.15f;

    #endregion

    EnergyCurve Curve { get { return mapGenTool.GetEnergyCurve(); } }
    List<CurveVertex> PolyPoints { get { return Curve.verts; } }
    List<Obstacle> Obstacles { get { return Curve.obstacles; } }
    Spline Spline { get { return mapGenTool.roadSpline; } }
    #endregion

    // Start is called before the first frame update
    void Start()
    {
        mapGenTool = new MapGenTool(
            object_road,
            object_bridgeAscent,
            object_bridge,
            object_bridgeDescent,
            object_crossing,
            object_ramp
        );  
    }

    void OnGUI()
    {
        // You can only call GUI functions in OnGUI()
        mapGenTool.style_Label = new GUIStyle(GUI.skin.label) { fontSize = 16 };
        mapGenTool.style_InlineLabel = new GUIStyle(mapGenTool.style_Label) { alignment = TextAnchor.MiddleLeft };
        mapGenTool.style_TextField = new GUIStyle(GUI.skin.textField) { fontSize = 16 };
        mapGenTool.style_Button = new GUIStyle(GUI.skin.button) { fontSize = 16 };
        mapGenTool.style_Slider = new GUIStyle(GUI.skin.horizontalSlider) { fontSize = 16 };
        mapGenTool.style_SliderThumb = new GUIStyle(GUI.skin.horizontalSliderThumb) { fontSize = 16 };
        mapGenTool.style_Toggle = new GUIStyle(GUI.skin.toggle) { fontSize = 16 };

        mapGenTool.UpdateUI();
    }

    void Update()
    {
        mapGenTool.MyUpdate();
    }

    #region Drawing
    private void OnPostRender()
    {
        if (!mat)
        {
            Debug.LogError("Please assign a Material");
            return;
        }

        GL.PushMatrix();
        mat.SetPass(0);

        #region Polyline
        if (PolyPoints != null)
        {

            List<Vector3> polyPoints3D = new();
            for (int i = 0; i < PolyPoints.Count; i++)
                polyPoints3D.Add(MyMath.xzToX0Z(PolyPoints[i].Position()));

            if (mapGenTool.showPolyPoints)
            {
                GL.Begin(GL.TRIANGLES);
                GL.Color(polyColor);
                for (int i = 0; i < polyPoints3D.Count; i++)
                {
                    DrawGLDisc(polyPoints3D[i], polylineDiscWidth);
                }
                GL.End();
            }

            if (mapGenTool.curveClosed)
                polyPoints3D.Add(polyPoints3D[0]);

            if (mapGenTool.showPolyLine)
            {
                GL.Begin(GL.QUADS);
                GL.Color(polyColor);
                for (int i = 0; i < polyPoints3D.Count - 1; i++)
                {
                    DrawGLLine(polyPoints3D[i], polyPoints3D[i + 1], polylineWidth);
                }
                GL.End();
            }
        }
        #endregion

        #region Obstacles
        if (mapGenTool.showObstacles && Obstacles != null)
        {

            foreach (Obstacle obs in Obstacles)
            {
                if (!obs.IsEnabled)
                    continue;

                List<Vector3> obstaclePoints3D = new();
                foreach (CurveVertex vert in obs.verts)
                    obstaclePoints3D.Add(MyMath.xzToX0Z(vert.Position()));

                if (mapGenTool.showPolyPoints)
                {
                    GL.Begin(GL.TRIANGLES);
                    GL.Color(obstacleColor);
                    for (int i = 0; i < obstaclePoints3D.Count; i++)
                    {
                        DrawGLDisc(obstaclePoints3D[i], obstacleDiscWidth);
                    }
                    GL.End();
                }

                if (mapGenTool.curveClosed)
                    obstaclePoints3D.Add(obstaclePoints3D[0]);

                if (mapGenTool.showPolyLine)
                {
                    GL.Begin(GL.QUADS);
                    GL.Color(obstacleColor);
                    for (int i = 0; i < obstaclePoints3D.Count - 1; i++)
                    {
                        DrawGLLine(obstaclePoints3D[i], obstaclePoints3D[i + 1], obstacleWidth);
                    }
                    GL.End();
                }
            }
        }
        #endregion

        #region Spline
        if (Spline != null && mapGenTool.showBezier)
        {
            for (int i = 0; i < Spline.NumSegments; i++)
            {
                try
                {
                    Spline.GetPointsInSegment(i);
                }
                catch (Exception)
                {
                    Spline.CalculateAllPointsInSegments();
                }

                Vector2[] segmentPoints = Spline.GetPointsInSegment(i);
                List<Vector3> segmentPoints3D = new();
                foreach (Vector2 splinePoint in segmentPoints)
                    segmentPoints3D.Add(MyMath.xzToX0Z(splinePoint));

                // lines to handles
                GL.Begin(GL.QUADS);
                GL.Color(lineColor);
                if (mapGenTool.showBezierHandles)
                {
                    DrawGLLine(segmentPoints3D[1], segmentPoints3D[0], bezierLineWidth);
                    DrawGLLine(segmentPoints3D[2], segmentPoints3D[3], bezierLineWidth);
                }
                GL.End();

                // bezier curve
                GL.Begin(GL.QUADS);
                GL.Color(splineColor);
                DrawBezierCurve(segmentPoints3D[0], segmentPoints3D[1], segmentPoints3D[2], segmentPoints3D[3], bezierWidth);
                GL.End();

                // handles
                GL.Begin(GL.TRIANGLES);
                GL.Color(anchorHandleColor);
                DrawGLDisc(segmentPoints3D[0], bezierHandleWidth);
                DrawGLDisc(segmentPoints3D[3], bezierHandleWidth);
                GL.Color(controlHandleColor);
                DrawGLDisc(segmentPoints3D[1], bezierHandleWidth);
                DrawGLDisc(segmentPoints3D[2], bezierHandleWidth);
                GL.End();
            }
        }
        #endregion

        GL.PopMatrix();
    }

    void DrawGLLine(Vector3 start, Vector3 end, float width)
    {
        Vector3 direction = (end - start).normalized;
        Vector3 cross = Vector3.Cross(direction, Vector3.up).normalized * width / 2.0f;

        GL.Vertex(start - cross);
        GL.Vertex(start + cross);
        GL.Vertex(end + cross);
        GL.Vertex(end - cross);
    }

    void DrawGLDisc(Vector3 center, float radius)
    {
        int resolution = 20;
        for (int i = 0; i < resolution; i++)
        {
            float angle = i * 2.0f * Mathf.PI / resolution;
            float nextAngle = (i + 1) * 2.0f * Mathf.PI / resolution;

            Vector3 vertex1 = center + new Vector3(Mathf.Cos(angle) * radius, 0, Mathf.Sin(angle) * radius);
            Vector3 vertex2 = center + new Vector3(Mathf.Cos(nextAngle) * radius, 0, Mathf.Sin(nextAngle) * radius);

            GL.Vertex(center);
            GL.Vertex(vertex2);
            GL.Vertex(vertex1);
        }
    }

    void DrawBezierCurve(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, float width, int resolution = 20)
    {
        Vector3 previousPoint = p0;

        for (int i = 1; i <= resolution; i++)
        {
            float t = i / (float)resolution;
            Vector3 currentPoint = Bezier.EvaluateCubic(p0, p1, p2, p3, t);

            DrawGLLine(previousPoint, currentPoint, width);
            previousPoint = currentPoint;
        }
    }
    #endregion


    private class MapGenTool : MapGenToolBase
    {
        private Vector2 scrollPosition;

        public float maxWidth = 400;
        public float contentHeight = 2000;

        public GUIStyle style_Label;
        public GUIStyle style_InlineLabel;
        public GUIStyle style_TextField;
        public GUIStyle style_Button;
        public GUIStyle style_Slider;
        public GUIStyle style_SliderThumb;
        public GUIStyle style_Toggle;

        Dictionary<int, bool> foldouts;

        EnergyCurve _curve;
        public override EnergyCurve GetEnergyCurve() => _curve ??= GeneratePolyline();
        protected override void SetEnergyCurve(EnergyCurve curve) => _curve = curve;

        readonly MeshTopologyEditorConfig _object_road;
        readonly MeshTopologyEditorConfig _object_bridgeAscent;
        readonly MeshTopologyEditorConfig _object_bridge;
        readonly MeshTopologyEditorConfig _object_bridgeDescent;
        readonly MeshTopologyEditorConfig _object_crossing;
        readonly MeshTopologyEditorConfig _object_ramp;
        protected override MeshTopologyEditorConfig GetObject_road() => _object_road;
        protected override MeshTopologyEditorConfig GetObject_bridgeAscent() => _object_bridgeAscent;
        protected override MeshTopologyEditorConfig GetObject_bridge() => _object_bridge;
        protected override MeshTopologyEditorConfig GetObject_bridgeDescent() => _object_bridgeDescent;
        protected override MeshTopologyEditorConfig GetObject_crossing() => _object_crossing;
        protected override MeshTopologyEditorConfig GetObject_ramp() => _object_ramp;

        public MapGenTool(
            MeshTopologyEditorConfig object_road,
            MeshTopologyEditorConfig object_bridgeAscent,
            MeshTopologyEditorConfig object_bridge,
            MeshTopologyEditorConfig object_bridgeDescent,
            MeshTopologyEditorConfig object_crossing,
            MeshTopologyEditorConfig object_ramp
        ) {
            _object_road = object_road;
            _object_bridgeAscent = object_bridgeAscent;
            _object_bridge = object_bridge;
            _object_bridgeDescent = object_bridgeDescent;
            _object_crossing = object_crossing;
            _object_ramp = object_ramp;

            foldouts = new();
        }

        public override void UpdateUI()
        {
            // Background
            GUI.Box(new Rect(0, 0, maxWidth, Screen.height), GUIContent.none);

            GUILayout.BeginHorizontal();
            GUILayout.Space(10);
            GUILayout.BeginVertical();
            GUILayout.Space(10);

            // Begin the ScrollView
            scrollPosition = GUI.BeginScrollView(
                new Rect(0, 0, maxWidth, Screen.height),
                scrollPosition,
                new Rect(0, 0, maxWidth - 20, contentHeight) // Subtracting 20 for the scrollbar width
            );

            base.UpdateUI();

            GUI.EndScrollView();
            GUILayout.EndVertical();
            GUILayout.EndHorizontal();
        }

        public override void RegenerateRoad()
        {
            Debug.Log("Not functional at the moment");
        }

        protected override void LoadMap()
        {
            Debug.Log("Not functional at the moment");
        }

        protected override void SaveMap()
        {
            Debug.Log("Not functional at the moment");
        }

        protected override void ExportPrefab()
        {
            Debug.Log("Not functional at the moment");
        }

        #region GUI methods
        protected override void CreateSection(string name, Action content)
        {
            CreateFoldout(name, content, 24);
        }

        protected override void CreateSubSection(string name, Action content)
        {
            CreateFoldout(name, content, 20);
        }

        protected override void CreateSubSubSection(string name, Action content)
        {
            CreateFoldout(name, content, 16);
        }

        protected override void CreateLabel(string name, int fontSize = 12)
        {
            GUILayout.Label(name, style_Label);
        }

        protected override void CreateTextField(string name, ref string input)
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label(name, style_InlineLabel);
            input = GUILayout.TextField(input, style_TextField);
            GUILayout.EndHorizontal();
        }

        protected override void CreateIntField(string name, ref int input)
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label(name, style_InlineLabel);
            input = int.TryParse(GUILayout.TextField(input.ToString(), style_TextField), out int result) ? result : input;
            GUILayout.EndHorizontal();
        }

        protected override void CreateFloatField(string name, ref float input)
        {
            GUILayout.BeginHorizontal();
            GUILayout.Label(name, style_InlineLabel);
            input = float.TryParse(GUILayout.TextField(input.ToString(), style_TextField), out float result) ? result : input;
            GUILayout.EndHorizontal();
        }

        protected override void CreateVector2Field(string name, ref Vector2 input)
        {
            CreateFloatField(name + " X", ref input.x);
            CreateFloatField(name + " Y", ref input.y);
        }

        protected override void CreateButton(string name, Action callback)
        {
            if (GUILayout.Button(name, style_Button)) callback();
        }

        protected override void CreateCheckbox(string name, ref bool trigger)
        {
            trigger = GUILayout.Toggle(trigger, name, style_Toggle);
        }

        protected override void CreateCheckbox_Dict<T>(string name, Dictionary<T, bool> dict, T key)
        {
            dict[key] = GUILayout.Toggle(dict[key], name, style_Toggle);
        }

        protected override void CreateSlider(ref int value, string text, int start, int end)
        {
            GUILayout.Label($"{text}: {value}");
            float rawValue = GUILayout.HorizontalSlider(value, start, end, style_Slider, style_SliderThumb);
            value = Mathf.RoundToInt(rawValue);
        }

        protected override void CreateSliderFloat(ref float value, string text, float start, float end, float snap = 0.5f)
        {
            GUILayout.Label($"{text}: {value:F1}");
            float rawValue = GUILayout.HorizontalSlider(value, start, end, style_Slider, style_SliderThumb);
            value = Mathf.Round(rawValue / snap) * snap;
        }

        protected override void CreateCurveField(string name, ref AnimationCurve curve)
        {
            throw new NotImplementedException();
        }

        protected override void CreateFoldout(string name, Action content, int fontSize, int bottomSpacing = 0)
        {
            GUIStyle foldoutStyle = new()
            {
                fontSize = fontSize,
                fontStyle = FontStyle.Bold,
                normal = {
                    textColor = Color.white
                },
                hover = {
                    background = Texture2D.whiteTexture,
                    textColor = Color.blue
                }
            };


            int foldoutHash = content.GetHashCode();

            foldouts[foldoutHash] = GUILayout.Toggle(
                foldouts.ContainsKey(foldoutHash) && foldouts[foldoutHash],
                name,
                foldoutStyle
            );

            if (foldouts[foldoutHash])
            {
                GUILayout.BeginHorizontal();
                GUILayout.Space(10);
                GUILayout.BeginVertical();

                content();

                GUILayout.EndVertical();
                GUILayout.EndHorizontal();
            }
        }

        protected override void CreateEnumSelection<T>(string name, T currentSelection, Action<T> onValueChanged)
        {
            // Create a section with the enum name
            GUILayout.Label(name, style_Label);

            // Get all values of the enum
            T[] enumValues = (T[])Enum.GetValues(typeof(T));

            // Create a button for each enum value
            foreach (T enumValue in enumValues)
            {
                // Determine if this is the currently selected value
                bool isSelected = currentSelection.Equals(enumValue);

                // Create a button with custom styling based on selection
                GUIStyle buttonStyle = new (style_Button)
                {
                    fontStyle = isSelected ? FontStyle.Bold : FontStyle.Normal,
                    normal = {
                        textColor = isSelected ? Color.green : Color.white
                    }
                };

                // Create the button
                if (GUILayout.Button(enumValue.ToString(), buttonStyle))
                {
                    onValueChanged(enumValue);
                }
            }
        }
        #endregion
    }
}
