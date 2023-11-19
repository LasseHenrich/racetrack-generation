using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Drawing : MonoBehaviour
{
    [SerializeField] Transform curveObjectsContainer;
    [SerializeField] Material mat;

    #region Colors

    // Note: Different colors don't work atm, probably because we only use one obstacle
    [SerializeField] Color controlHandleColor = new(0.2f, 0.5f, 0.2f);
    [SerializeField] Color anchorHandleColor = new(255f / 255f, 225f / 255f, 67f / 255f);
    [SerializeField] Color lineColor = new(0, 0, 0);
    [SerializeField] Color splineColor = new(0.2f, 0.2f, 0.5f);
    [SerializeField] Color polyColor = new(0.4f, 0.6f, 1f);
    [SerializeField] Color obstacleColor = new(1f, 0f, 0f);

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

    ControlWindow ctrlWindow;
    EnergyCurve Curve { get { return ctrlWindow.Curve; } }
    List<CurveVertex> PolyPoints { get { return Curve.verts; } }
    List<Obstacle> Obstacles { get { return Curve.obstacles; } }

    void Start()
    {
        ctrlWindow = FindObjectOfType<ControlWindow>();
    }

    private void OnPostRender()
    {
        if (!mat)
        {
            Debug.LogError("Please assign a Material");
            return;
        }

        GL.PushMatrix();
        if (!mat.SetPass(0))
            return;

        #region Polyline
        if (PolyPoints != null)
        {
            GL.Color(polyColor);

            List<Vector3> polyPoints = new();
            for (int i = 0; i < PolyPoints.Count; i++)
                polyPoints.Add(PolyPoints[i].Position());
            if (ctrlWindow.curveClosed) polyPoints.Add(PolyPoints[0].Position());

            if (ctrlWindow.showPolyLine)
                DrawGLLinesFromList(polyPoints);

            if (ctrlWindow.showPolyPoints)
                DrawGLDiscsFromList(polyPoints, polylineWidth);
        }
        #endregion

        #region Obstacles
        if (Obstacles != null && ctrlWindow.showObstacles)
        {
            GL.Color(obstacleColor);

            foreach (Obstacle obs in Obstacles)
            {
                List<Vector3> obstaclePoints = new();
                foreach (CurveVertex vert in obs.verts)
                    obstaclePoints.Add(vert.Position());
                obstaclePoints.Add(obs.verts[0].Position());

                if (ctrlWindow.showPolyLine)
                    DrawGLLinesFromList(obstaclePoints);

                if (ctrlWindow.showPolyPoints)
                    DrawGLDiscsFromList(obstaclePoints, polylineWidth);
            }

        }
        #endregion

        GL.PopMatrix();
    }

    void DrawGLLinesFromList(List<Vector3> list)
    {
        GL.Begin(GL.QUADS);
        for (int i = 0; i < list.Count - 1; i++)
            DrawGLLine(list[i], list[i + 1], polylineWidth);
        GL.End();
    }

    void DrawGLDiscsFromList(List<Vector3> list, float radius)
    {
        GL.Begin(GL.TRIANGLES);
        for (int i = 0; i < list.Count - 1; i++)
            DrawGLDisc(list[i], radius);
        GL.End();
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
}
