using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class PlaneObstacle : Obstacle
{
    public Vector3 pointA, pointB, pointC;
    private Vector3 planeNormal;
    private float p_exp;

    public PlaneObstacle(Vector3 pointA, Vector3 pointB, Vector3 pointC, float p_exp) : base()
    {
        this.pointA = pointA;
        this.pointB = pointB;
        this.pointC = pointC;
        this.p_exp = p_exp;

        // Calculate the plane normal using the cross product of two edge vectors
        planeNormal = Vector3.Cross(pointB - pointA, pointC - pointA).normalized;
    }

    public override void AddGradient(EnergyCurve curve, Matrix<float> gradient)
    {
        int numVerts = curve.NumVerts();

        for (int i = 0; i < numVerts; i++)
        {
            CurveVertex v = curve.verts[i];
            Vector3 toPoint = VectorToClosestPoint(v.Position());
            if (toPoint == Vector3.zero)
                continue;

            float dist = toPoint.magnitude;
            Vector3 grad = toPoint.normalized * p_exp / Mathf.Pow(dist, p_exp + 1);

            CurveGenUtils.AddToRow(gradient, v.GlobalIndex(), grad);
        }
    }

    public override float ComputeEnergy(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        float sumE = 0;

        for (int i = 0; i < numVerts; i++)
        {
            Vector3 toPoint = VectorToClosestPoint(curve.verts[i].Position());
            if (toPoint == Vector3.zero)
                continue;

            float dist = toPoint.magnitude;
            sumE += 1f / Mathf.Pow(dist, p_exp);
        }

        return sumE;
    }

    private Vector3 ClosestPointOnPlane(Vector3 point)
    {
        // Calculate the projection of the point onto the plane
        float d = Vector3.Dot(planeNormal, pointA - point);
        return point + planeNormal * d;
    }

    private Vector3 VectorToClosestPoint(Vector3 point)
    {
        Vector3 nearest = ClosestPointOnPlane(point);
        return nearest - point;
    }
}
