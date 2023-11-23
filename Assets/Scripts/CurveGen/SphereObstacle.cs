using MathNet.Numerics.LinearAlgebra;
using UnityEngine;

public class SphereObstacle : Obstacle
{
    public Vector3 center;
    public float radius, p_exp;

    public SphereObstacle(Vector3 center, float radius, float p_exp) : base()
    {
        this.center = center;
        this.radius = radius;
        this.p_exp = p_exp;
    }

    public override void AddGradient(EnergyCurve curve, Matrix<float> gradient)
    {
        int numVerts = curve.NumVerts();

        for (int i = 0; i < numVerts; i++)
        {
            CurveVertex v = curve.verts[i];

            Vector3 toPoint = VectorToClosestPoint(v);
            if (toPoint == Vector3.zero)
                continue;

            float dist = toPoint.magnitude;
            toPoint /= dist;
            Vector3 grad = toPoint * p_exp / Mathf.Pow(dist, p_exp + 1);

            CurveGenUtils.AddToRow(gradient, v.GlobalIndex(), grad);
        }
    }

    public override float ComputeEnergy(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        float sumE = 0;

        for (int i = 0; i < numVerts; i++)
        {
            Vector3 toPoint = VectorToClosestPoint(curve.verts[i]);
            if (toPoint == Vector3.zero)
                continue;

            float dist = toPoint.magnitude;
            sumE += 1f / Mathf.Pow(dist, p_exp); // todo: refactor two lines to use sqrMagnitude
        }

        return sumE;
    }

    private Vector3 ClosestPoint(Vector3 input)
    {
        Vector3 dir = (input - center).normalized;
        return center + radius * dir;
    }

    private Vector3 VectorToClosestPoint(CurveVertex vert)
    {
        Vector3 pos_i = vert.Position();
        // If we're very close to the center of the sphere, gradient is 0
        if ((pos_i - center).sqrMagnitude < 1e-12)
            return Vector3.zero;
        // Find the closest point on the plane
        Vector3 nearest = ClosestPoint(pos_i);
        // Simulate an energy contribution of 1 / r^(b - a)
        return nearest - pos_i;
    }
}
