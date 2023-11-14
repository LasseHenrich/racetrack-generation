using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using UnityEngine;

public class Obstacle : Curve // ToDo: Make inherit from Potential
{
    public readonly float p_exp, weight, radius;
    public readonly int numPoints;
    readonly BVHNode2D root;
    bool enabled;
    public Vector2 center; // For Serialization

    public bool IsEnabled { get { return enabled; } }

    public Obstacle(float p_exp, float weight, int numPoints, float radius, Vector2 center)
    {
        this.p_exp = p_exp;
        this.weight = weight;
        this.numPoints = numPoints;
        this.radius = radius;
        this.center = center;
        enabled = true;

        curveClosed = true;

        List<Vector2> vertPositions = new List<Vector2>();
        for (int i = 0; i < numPoints; i++)
        {
            float alpha = (i / (float)numPoints) * Mathf.PI * 2f;
            vertPositions.Add(center + radius * new Vector2(Mathf.Cos(alpha), Mathf.Sin(alpha)));
        }
        InitVertsEdgesFromPositions(vertPositions);

        root = CreateBVHFromVerts();
        //Debug.Log("num childs: " + bvh.TotalLeafCount());
    }

    private BVHNode2D CreateBVHFromVerts()
    {
        int numVerts = verts.Count;
        List<VertexBody4D> vertBodies = new();

        for (int i = 0; i < numVerts; i++)
        {
            VertexBody4D currBody = VertToBody(verts[i]);
            vertBodies.Add(currBody);
        }

        BVHNode2D tree = new(vertBodies, 0, null, false);
        tree.RecomputeCentersOfMass(this);
        BVHNode2D.globalID = 0;
        tree.RecursivelyAssignIDs();
        return tree;
    }

    private VertexBody4D VertToBody(CurveVertex v)
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
            Vector2 pos = curve.verts[i].Position();
            Vector2 force = AccumulateForce(root, pos);
            CurveGenUtils.AddToRow(gradient, i, force * weight);
        }
    }

    private Vector2 AccumulateForce(BVHNode2D node, Vector2 point)
    {
        if (node.isEmpty)
            return Vector2.zero;

        if (node.isLeaf)
            return BodyForce(node, point);

        if (node.ShouldUseCell(point))
            return BodyForce(node, point);

        Vector2 total = Vector2.zero;
        foreach (BVHNode2D child in node.children)
            total += AccumulateForce(child, point);
        return total;
    }

    private Vector2 BodyForce(BVHNode2D node, Vector2 point)
    {
        Vector2 center = node.centerOfMass;
        float mass = node.totalMass;

        Vector2 toPoint = center - point;
        float dist = toPoint.magnitude;
        toPoint /= dist;

        Vector2 grad_i = toPoint * p_exp / Mathf.Pow(dist, p_exp + 1);

        return mass * grad_i;
    }

    internal float ComputeEnergy(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        float sumE = 0;
        for (int i = 0; i < numVerts; i++)
        {
            Vector2 pos = curve.verts[i].Position();
            sumE += AccumulateEnergy(root, pos);
        }
        return weight * sumE;
    }

    private float AccumulateEnergy(BVHNode2D node, Vector2 point)
    {
        if (node.isEmpty)
            return 0;

        if (node.isLeaf)
            return BodyEnergy(node, point);

        if (node.ShouldUseCell(point))
            return BodyEnergy(node, point);

        float total = 0;
        foreach (BVHNode2D child in node.children)
            total += AccumulateEnergy(child, point);
        return total;
    }

    private float BodyEnergy(BVHNode2D node, Vector2 point)
    {
        Vector2 center = node.centerOfMass;
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
