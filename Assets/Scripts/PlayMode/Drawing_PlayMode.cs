using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Drawing_PlayMode : MonoBehaviour
{
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

    ControlWindow_PlayMode ctrlWindow;
    EnergyCurve Curve { get { return ctrlWindow.Curve; } }
    List<CurveVertex> PolyPoints { get { return Curve.verts; } }
    List<Obstacle> Obstacles { get { return Curve.obstacles; } }
    Spline Spline { get { return ctrlWindow.roadSpline; } }

    void Start()
    {
        ctrlWindow = FindObjectOfType<ControlWindow_PlayMode>();
    }

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

            if (ctrlWindow.showPolyPoints)
            {
                GL.Begin(GL.TRIANGLES);
                GL.Color(polyColor);
                for (int i = 0; i < polyPoints3D.Count; i++)
                {
                    DrawGLDisc(polyPoints3D[i], polylineDiscWidth);
                }
                GL.End();
            }

            if (ctrlWindow.curveClosed)
                polyPoints3D.Add(polyPoints3D[0]);

            if (ctrlWindow.showPolyLine)
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
        if (ctrlWindow.showObstacles && Obstacles != null)
        {

            foreach (Obstacle obs in Obstacles)
            {
                if (!obs.IsEnabled)
                    continue;

                List<Vector3> obstaclePoints3D = new();
                foreach (CurveVertex vert in obs.verts)
                    obstaclePoints3D.Add(MyMath.xzToX0Z(vert.Position()));

                if (ctrlWindow.showPolyPoints)
                {
                    GL.Begin(GL.TRIANGLES);
                    GL.Color(obstacleColor);
                    for (int i = 0; i < obstaclePoints3D.Count; i++)
                    {
                        DrawGLDisc(obstaclePoints3D[i], obstacleDiscWidth);
                    }
                    GL.End();
                }

                if (ctrlWindow.curveClosed)
                    obstaclePoints3D.Add(obstaclePoints3D[0]);

                if (ctrlWindow.showPolyLine)
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
        if (Spline != null && ctrlWindow.showBezier)
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
                if (ctrlWindow.showBezierHandles)
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

}
