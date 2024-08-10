using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Drawing_PlayMode : MonoBehaviour
{
    [SerializeField] Material mat;

    #region Colors

    [SerializeField] Color controlHandleColor = new(0.2f, 0.5f, 0.2f);
    [SerializeField] Color anchorHandleColor = new(255f / 255f, 225f / 255f, 67f / 255f);
    [SerializeField] Color lineColor = new(0, 0, 0);
    [SerializeField] Color splineColor = new(0.2f, 0.2f, 0.5f);
    [SerializeField] Color polyColor = new(0.4f, 0.6f, 1f);

    #endregion

    #region Sizes

    //[SerializeField] float bezierWidth = 8f; // 20f
    //[SerializeField] float bezierLineWidth = 3f; // 5f
    //[SerializeField] float bezierHandleWidth = 1.4f;//1.0f;//0.6f;
    [SerializeField] float polylineWidth = 4f;
    [SerializeField] float polylineDiscWidth = 0.15f;
    //[SerializeField] float obstacleWidth = 2f;
    //[SerializeField] float obstacleDiscWidth = 0.15f;

    #endregion

    ControlWindow_PlayMode ctrlWindow;
    EnergyCurve Curve { get { return ctrlWindow.Curve; } }
    List<CurveVertex> PolyPoints { get { return Curve.verts; } }

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
            GL.Begin(GL.QUADS);

            GL.Color(polyColor);
            List<Vector3> polyPoints3D = new();
            for (int i = 0; i < PolyPoints.Count; i++)
                polyPoints3D.Add(MyMath.xzToX0Z(PolyPoints[i].Position()));
            if (ctrlWindow.curveClosed) polyPoints3D.Add(MyMath.xzToX0Z(PolyPoints[0].Position()));

            if (ctrlWindow.showPolyLine)
            {
                for (int i = 0; i < polyPoints3D.Count - 1; i++)
                {
                    DrawGLLine(polyPoints3D[i], polyPoints3D[i + 1], polylineWidth);
                }
            }

            GL.End();

            if (ctrlWindow.showPolyPoints)
            {
                GL.Begin(GL.TRIANGLES);
                for (int i = 0; i < PolyPoints.Count; i++)
                {
                    DrawGLDisc(MyMath.xzToX0Z(PolyPoints[i].Position()), polylineDiscWidth);
                }
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
}
