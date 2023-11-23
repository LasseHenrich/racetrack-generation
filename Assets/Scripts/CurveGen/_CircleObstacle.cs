using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using UnityEngine;

public class _CircleObstacle : Curve // ToDo: Make inherit from Potential
{
    public readonly float p_exp, weight, radius;
    public readonly int numPoints;
    readonly BVHNode3D root;
    bool enabled;
    public Vector3 center; // For Serialization

    public bool IsEnabled { get { return enabled; } }

    public _CircleObstacle(float p_exp, float weight, int numPoints, float radius, Vector3 center)
    {
        this.p_exp = p_exp;
        this.weight = weight;
        this.numPoints = numPoints;
        this.radius = radius;
        this.center = center;
        enabled = true;

        curveClosed = true;

        List<Vector3> vertPositions = new();
        for (int i = 0; i < numPoints; i++)
        {
            float alpha = (i / (float)numPoints) * Mathf.PI * 2f;
            vertPositions.Add(center + radius * new Vector3(Mathf.Cos(alpha), 0, Mathf.Sin(alpha)));
        }
        InitVertsEdgesFromPositions(vertPositions);

        root = CreateBVHFromVerts();
        //Debug.Log("num childs: " + bvh.TotalLeafCount());
    }

    private BVHNode3D CreateBVHFromVerts()
    {
        int numVerts = verts.Count;
        List<VertexBody6D> vertBodies = new();

        for (int i = 0; i < numVerts; i++)
        {
            VertexBody6D currBody = VertToBody(verts[i]);
            vertBodies.Add(currBody);
        }

        BVHNode3D tree = new(vertBodies, 0, null, false);
        tree.RecomputeCentersOfMass(this);
        BVHNode3D.globalID = 0;
        tree.RecursivelyAssignIDs();
        return tree;
    }

    private VertexBody6D VertToBody(CurveVertex v)
    {
        PosTan pt = new(v.Position(), v.Tangent());
        float mass = v.AvgLength();

        return new(pt, mass, v.GlobalIndex(), BodyType.Vertex);
    }

    public void AddGradient(EnergyCurve curve, Matrix<float> gradient)
    {
        int numVerts = curve.NumVerts();
        for (int i = 0; i < numVerts; i++)
        {
            Vector3 pos = curve.verts[i].Position();
            Vector3 force = AccumulateForce(root, pos);
            CurveGenUtils.AddToRow(gradient, i, force * weight);
        }
    }

    private Vector3 AccumulateForce(BVHNode3D node, Vector3 point)
    {
        if (node.isEmpty)
            return Vector3.zero;

        if (node.isLeaf)
            return BodyForce(node, point);

        if (node.ShouldUseCell(point))
            return BodyForce(node, point);

        Vector3 total = Vector3.zero;
        foreach (BVHNode3D child in node.children)
            total += AccumulateForce(child, point);
        return total;
    }

    private Vector3 BodyForce(BVHNode3D node, Vector3 point)
    {
        Vector3 center = node.centerOfMass;
        float mass = node.totalMass;

        Vector3 toPoint = center - point;
        float dist = toPoint.magnitude;
        toPoint /= dist;

        Vector3 grad_i = toPoint * p_exp / Mathf.Pow(dist, p_exp + 1);

        return mass * grad_i;
    }

    internal float ComputeEnergy(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        float sumE = 0;
        for (int i = 0; i < numVerts; i++)
        {
            Vector3 pos = curve.verts[i].Position();
            sumE += AccumulateEnergy(root, pos);
        }
        return weight * sumE;
    }

    private float AccumulateEnergy(BVHNode3D node, Vector3 point)
    {
        if (node.isEmpty)
            return 0;

        if (node.isLeaf)
            return BodyEnergy(node, point);

        if (node.ShouldUseCell(point))
            return BodyEnergy(node, point);

        float total = 0;
        foreach (BVHNode3D child in node.children)
            total += AccumulateEnergy(child, point);
        return total;
    }

    private float BodyEnergy(BVHNode3D node, Vector3 point)
    {
        Vector3 center = node.centerOfMass;
        float mass = node.totalMass;
        float distance = (center - point).magnitude;
        return 1.0f / Mathf.Pow(distance, p_exp);
    }

    public void Disable()
    {
        enabled = false;
    }

    public void Enable()
    {
        enabled = true;
    }
}
