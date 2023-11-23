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
                if (obs is SphereObstacle sphereObstacle)
                {
                    DrawGLSphere(sphereObstacle.center, sphereObstacle.radius);
                }
                else if (obs is PlaneObstacle planeObstacle)
                {
                    DrawGLPlane(planeObstacle.pointA, planeObstacle.pointB, planeObstacle.pointC);
                }
                else
                {
                    Debug.LogWarning("Drawing not yet implemented for obstacle type");
                }
            }

        }
        #endregion

        GL.PopMatrix();
    }

    void DrawGLLinesFromList(List<Vector3> list)
    {
        for (int i = 0; i < list.Count - 1; i++)
            DrawGLLine(list[i], list[i + 1], polylineWidth);
    }

    void DrawGLDiscsFromList(List<Vector3> list, float radius)
    {
        for (int i = 0; i < list.Count - 1; i++)
            DrawGLDisc(list[i], radius);
    }

    void DrawGLLine(Vector3 start, Vector3 end, float width)
    {
        GL.Begin(GL.QUADS);

        Vector3 direction = (end - start).normalized;
        Vector3 cross = Vector3.Cross(direction, Vector3.up).normalized * width / 2.0f;

        GL.Vertex(start - cross);
        GL.Vertex(start + cross);
        GL.Vertex(end + cross);
        GL.Vertex(end - cross);

        GL.End();
    }

    void DrawGLDisc(Vector3 center, float radius)
    {
        GL.Begin(GL.TRIANGLES);

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

        GL.End();
    }

    void DrawGLSphere(Vector3 center, float radius)
    {
        GL.Begin(GL.TRIANGLES);

        int resolution = 10;
        for (int i = 0; i <= resolution; i++)
        {
            float angle = i * Mathf.PI / resolution;
            Vector3 discCenter = new(center.x, Mathf.Sin(angle) * radius + center.y, center.z);
            float discRadius = Mathf.Cos(angle) * radius;

            DrawGLDisc(discCenter, discRadius);
        }

        GL.End();
    }

    void DrawGLPlane(Vector3 A, Vector3 B, Vector3 C)
    {
        GL.Begin(GL.TRIANGLES);

        Vector3 midpoint = (A + C) / 2;
        Vector3 vectorToAdjacent = B - midpoint;
        List<Vector3> possibleDs = new()
        {
            A + B - C,
            A - B + C,
            -A + B + C
        };

        Vector3 D = possibleDs.Find((possD) => possD.sqrMagnitude <= A.sqrMagnitude && possD.sqrMagnitude <= B.sqrMagnitude && possD.sqrMagnitude <= C.sqrMagnitude);

        Debug.Log($"{A}, {B}, {C}, {D}");

        // Draw the two triangles that make up the plane
        GL.Vertex(A);
        GL.Vertex(B);
        GL.Vertex(C);

        GL.Vertex(B);
        GL.Vertex(D);
        GL.Vertex(C);

        GL.End();
    }
}
