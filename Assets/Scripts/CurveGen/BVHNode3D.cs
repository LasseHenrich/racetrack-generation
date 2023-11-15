using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BVHNode3D
{
    public static int globalID;
    public int thisNodeID;
    public int numNodes;
    public float totalMass;
    public Vector3 centerOfMass;

    private Vector3 averageTangent;
    private readonly List<int> clusterIndices = new();
    private readonly BVHNode3D bvhRoot;
    private readonly Vector<float> fullMasses;

    // Fields for use by matrix-vector products.
    // Not used by any BVH functions.

    public List<BVHNode3D> children = new();
    private readonly VertexBody6D body = new();

    int numElements;

    readonly float splitPoint;
    public readonly bool isEmpty;
    public readonly bool isLeaf;
    private PosTan minCoords;
    private PosTan maxCoords;
    private const float thresholdTheta = 0.25f;

    public BVHNode3D(List<VertexBody6D> points, int axis, BVHNode3D root, bool splitTangents)
    {
        if (root == null) bvhRoot = this;
        else bvhRoot = root;

        if (points.Count == 0)
        {
            isLeaf = false;
            isEmpty = true;
            totalMass = 0;
            numElements = 0;
            // This is just to find bugs more easily
            body.elementIndex = -999;
        }
        // If there's only one point, then we don't need to do anything
        // except set the fields from that one point
        else if (points.Count == 1)
        {
            body = points[0];
            isLeaf = true;
            isEmpty = false;
            totalMass = body.mass;
            centerOfMass = body.pt.position;
            averageTangent = body.pt.tangent;
            minCoords = body.pt;
            maxCoords = body.pt;
            numElements = 1;

            clusterIndices.Add(body.elementIndex);
        }
        else
        {
            // Reserve space for splitting the points into lesser and greater
            int numPoints = points.Count;
            List<VertexBody6D> lesserPoints = new();
            List<VertexBody6D> greaterPoints = new();

            // Compute the plane over which to split the points;
            splitPoint = AxisSplittingPlane(points, axis);

            // split the points over the median
            for (int i = 0; i < numPoints; i++)
            {
                float coord = GetCoordFromBody(points[i], axis);

                if (coord <= splitPoint)
                    lesserPoints.Add(points[i]);
                else
                    greaterPoints.Add(points[i]);
            }

            // Compute the middle extents for the two halves

            // Recursively construct children
            int nextAxis;
            if (splitTangents)
                nextAxis = NextAxis(axis);
            else
                nextAxis = NextSpatialAxis(axis);

            BVHNode3D nextRoot = root ?? null;
            BVHNode3D lesserNode = new(lesserPoints, nextAxis, nextRoot, splitTangents);
            BVHNode3D greaterNode = new(greaterPoints, nextAxis, nextRoot, splitTangents);

            children.Add(lesserNode);
            children.Add(greaterNode);

            clusterIndices.AddRange(lesserNode.clusterIndices);
            clusterIndices.AddRange(greaterNode.clusterIndices);

            if (root == null)
            {
                fullMasses = Vector<float>.Build.Dense(points.Count);
                for (int i = 0; i < points.Count; i++)
                {
                    fullMasses[i] = points[i].mass;
                }
            }

            isLeaf = false;
            isEmpty = false;
        }
    }

    private int NextSpatialAxis(int axis)
    {
        if (axis < 3) return (axis + 1) % 3;
        else return 0;
    }

    private int NextAxis(int axis)
    {
        return (axis + 1) % 6;
    }

    private float AxisSplittingPlane(List<VertexBody6D> points, int axis)
    {
        int numPoints = points.Count;
        List<float> coords = new List<float>();

        for (int i = 0; i < numPoints; i++)
            coords.Add(GetCoordFromBody(points[i], axis));

        Sorting.Sort(coords); // Sort in non-descending order

        int splitIndex = -1;
        float minWidths = float.MaxValue;

        for (int i = 0; i < numPoints; i++)
        {
            float width1 = coords[i] - coords[0];
            float width2 = (i == numPoints - 1) ? 0 : coords[numPoints - 1] - coords[i + 1];

            float sumSquares = width1 * width1 + width2 * width2;
            if (sumSquares < minWidths)
            {
                minWidths = sumSquares;
                splitIndex = i;
            }
        }

        float splitPoint = (coords[splitIndex] + coords[splitIndex + 1]) * 0.5f;
        return splitPoint;
    }

    private float GetCoordFromBody(VertexBody6D body, int axis)
    {
        return axis switch
        {
            0 => body.pt.position.x,
            1 => body.pt.position.y,
            2 => body.pt.position.z,
            3 => body.pt.tangent.x,
            4 => body.pt.tangent.y,
            5 => body.pt.tangent.z,
            _ => throw new ArgumentException("Invalid axis passed to GetCoordFromBody"),
        };
    }

    internal void RecomputeCentersOfMass(Curve network)
    {
        if (isEmpty)
        {
            totalMass = 0;
            numElements = 0;
        }
        // For a leaf, just set centers and bounds from the one body
        else if (isLeaf)
        {
            SetLeafData(network);
            numElements = 1;
        }
        else
        {
            // Recursively compute bounds for all children
            for (int i = 0; i < children.Count; i++)
                children[i].RecomputeCentersOfMass(network);

            minCoords = children[0].minCoords;
            maxCoords = children[0].maxCoords;

            totalMass = 0;
            centerOfMass = Vector3.zero;
            averageTangent = Vector3.zero;

            // Accumulate max / min over all nonempty children
            for (int i = 0; i < children.Count; i++)
            {
                if (!children[i].isEmpty)
                {
                    minCoords = PosTanMin(children[i].minCoords, minCoords);
                    maxCoords = PosTanMax(children[i].maxCoords, maxCoords);

                    totalMass += children[i].totalMass;
                    centerOfMass += children[i].centerOfMass * children[i].totalMass;
                    averageTangent += children[i].averageTangent * children[i].totalMass;
                }
            }

            centerOfMass /= totalMass;
            averageTangent /= totalMass;

            averageTangent = averageTangent.normalized;

            numElements = 0;
            for (int i = 0; i < children.Count; i++)
                numElements += children[i].numElements;
        }
    }

    private PosTan PosTanMin(PosTan v1, PosTan v2)
    {
        return new PosTan(CurveGenUtils.VectorMin(v1.position, v2.position), CurveGenUtils.VectorMin(v1.tangent, v2.tangent));
    }

    private PosTan PosTanMax(PosTan v1, PosTan v2)
    {
        return new PosTan(CurveGenUtils.VectorMax(v1.position, v2.position), CurveGenUtils.VectorMax(v1.tangent, v2.tangent));
    }

    //private void SetLeafData<T>(T network)
    //{
    //    throw new ArgumentException("Called SetLeafData with type " + network.GetType() + ", which is not supported");
    //}

    private void SetLeafData(Curve curve)
    {
        if (curve == null)
            throw new ArgumentNullException("called SetLeafData with network == null");

        if (body.type == BodyType.Vertex)
        {
            CurveVertex p = curve.verts[body.elementIndex];
            body.mass = p.AvgLength();
            body.pt.position = p.Position();
            body.pt.tangent = p.Tangent();
        }
        else if (body.type == BodyType.Edge)
        {
            CurveEdge p1 = curve.edges[body.elementIndex];

            // Mass of an edge is its length
            body.mass = p1.Length();
            // Use midpoint as center of mass
            body.pt.position = p1.Midpoint();
            // Tangent direction is normalized edge vector
            body.pt.tangent = p1.Tangent();
        }

        totalMass = body.mass;
        centerOfMass = body.pt.position;
        averageTangent = body.pt.tangent;

        minCoords = new PosTan(body.pt.position, body.pt.tangent);
        maxCoords = minCoords;
    }

    internal void RecursivelyAssignIDs()
    {
        thisNodeID = globalID++;
        foreach (BVHNode3D child in children)
            child.RecursivelyAssignIDs();
    }

    internal void AccumulateTPEGradient(Matrix<float> gradients, CurveVertex i_pt, EnergyCurve curve, float alpha, float beta)
    {
        if (isEmpty)
            return;
        else if (isLeaf)
        {
            // With a vertex, we add gradient terms the same way as usual
            if (body.type == BodyType.Vertex)
            {
                // If this is a leaf, then it only has one vertex in it, so just use it
                CurveVertex j_pt = curve.verts[body.elementIndex];
                // Don't sum if it's the same vertex
                if (j_pt == i_pt) return;

                // Add i and neighbors of i
                List<CurveVertex> i_pts = new();
                i_pts.Add(i_pt);
                for (int e = 0; e < i_pt.NumEdges(); e++)
                    i_pts.Add(i_pt.Edge(e).Opposite(i_pt));

                // Add j and neighbors of j
                List<CurveVertex> j_pts = new();
                j_pts.Add(j_pt);
                for (int e = 0; e < j_pt.NumEdges(); e++)
                    j_pts.Add(j_pt.Edge(e).Opposite(j_pt));

                foreach (CurveVertex i_n in i_pts)
                {
                    CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(i_pt, j_pt, i_n));
                    bool noOverlap = true;
                    // Avoid double-counting on the reverse termns with these checks
                    foreach (CurveVertex j_n in j_pts)
                    {
                        if (i_n.Equals(j_n)) noOverlap = false;
                    }
                    if (noOverlap)
                        CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(j_pt, i_pt, i_n));
                }
            }
            // With an edge, we have to pass in a TangentMassPoint instead
            else if (body.type == BodyType.Edge)
            {
                Vector3 tangent = body.pt.tangent;
                tangent = tangent.normalized;
                CurveEdge edge = curve.edges[body.elementIndex];
                CurveVertex j1 = edge.GetPrevVertex();
                CurveVertex j2 = edge.GetNextVertex();

                // Add i and neighbors of i
                List<CurveVertex> i_pts = new();
                i_pts.Add(i_pt);
                for (int e = 0; e < i_pt.NumEdges(); e++)
                    i_pts.Add(i_pt.Edge(e).Opposite(i_pt));

                TangentMassPoint jm = new(tangent, body.mass, body.pt.position, j1, j2);

                if (i_pt != j1 && i_pt != j2)
                {
                    foreach (CurveVertex i_n in i_pts)
                    {
                        CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(i_pt, jm, i_n));
                        CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(jm, i_pt, i_n));
                    }
                }
            }
        }
        else
        {
            if (ShouldUseCell(i_pt.Position()))
            {
                Vector3 tangent = averageTangent;
                tangent = tangent.normalized;
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j = new(tangent, totalMass, centerOfMass, null, null);

                // Add i and neighbors of i
                List<CurveVertex> i_pts = new();
                i_pts.Add(i_pt);
                for (int e = 0; e < i_pt.NumEdges(); e++)
                    i_pts.Add(i_pt.Edge(e).Opposite(i_pt));

                // Differentiate both terms for previous, middle and next
                foreach (CurveVertex i_n in i_pts)
                {
                    CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(i_pt, j, i_n));
                    CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TPE.GetInstance.TpeGrad(j, i_pt, i_n));
                }
            }
            else
            {
                // Otherwise we continue recursively traversing the tree
                for (int i = 0; i < children.Count; i++)
                {
                    if (children[i] != null)
                        children[i].AccumulateTPEGradient(gradients, i_pt, curve, alpha, beta);
                }
            }
        }
    }

    public bool ShouldUseCell(Vector3 vertPos)
    {
        // ToDo: Take into account some tangent-related criteria?
        float f = (centerOfMass - vertPos).magnitude;
        bool shouldUseCell = NodeRatio(f) < thresholdTheta;
        //if (shouldUseCell) Debug.Log("Summarizing " + TotalLeafCount() + " leafs");
        return shouldUseCell;
    }

    private float NodeRatio(float f)
    {
        // Compute diagonal distance from corner to corner
        Vector3 diag = maxCoords.position - minCoords.position;
        return diag.magnitude / f;
    }

    internal float TotalEnergy(EnergyCurve curve, BVHNode3D root)
    {
        int numVerts = curve.NumVerts();
        float fullSum = 0;

        // Loop over all vertices and add up energy contributions
        for (int i = 0; i < numVerts; i++)
        {
            CurveVertex i_pt = curve.verts[i];
            fullSum += root.AccumulateVertexEnergy(i_pt, curve);
        }
        return fullSum;
    }

    private float AccumulateVertexEnergy(CurveVertex i_pt, EnergyCurve curve)
    {
        if (isEmpty)
            return 0;

        if (isLeaf)
        {
            // If this is a leaf, then it only has one element in it, so just use it
            if (body.type == BodyType.Vertex)
            {
                // If the element is a vertex, just get that vertex
                CurveVertex j_pt = curve.verts[body.elementIndex];
                // Don't sum if it's the same vertex
                if (j_pt.Equals(i_pt)) return 0;
                // Add contribution to energy at i from vertex j
                return TPE.GetInstance.TpePair(i_pt, j_pt);
            }
            else if (body.type == BodyType.Edge)
                // Otherwise the element is an edge, so we use its midpoint (stored in the body)
                return TPE.GetInstance.TpePairPts(i_pt.Position(), body.pt.position, i_pt.Tangent(), i_pt.AvgLength(), body.mass);
        }
        else
        {
            if (ShouldUseCell(i_pt.Position()))
                // This cell is far enough away that we can treat it as a single body
                return BodyEnergyEvaluation(i_pt);

            // Otherise, we continue recursively traversing the tree
            float result = 0;
            for (int i = 0; i < children.Count; i++)
            {
                if (children[i] != null)
                    result += children[i].AccumulateVertexEnergy(i_pt, curve);
            }
            return result;
        }
        return 0;
    }

    private float BodyEnergyEvaluation(CurveVertex i_pt)
    {
        Vector3 tangent = averageTangent;
        tangent = tangent.normalized;
        return TPE.GetInstance.TpePairPts(i_pt.Position(), centerOfMass, tangent, i_pt.AvgLength(), totalMass);
    }

    public int TotalLeafCount()
    {
        if (isLeaf) return 1;

        int childCount = 0;
        foreach (BVHNode3D child in children)
            childCount += child.TotalLeafCount();

        return childCount;
    }
}

public enum TPEPointType
{
    Cluster,
    Point,
    Edge
}

public class TangentMassPoint
{
    public Vector3 tangent;
    public float mass;
    public Vector3 point;
    public CurveVertex curvePt;
    public CurveVertex curvePt2;

    public TangentMassPoint(Vector3 tangent, float mass, Vector3 point, CurveVertex curvePt, CurveVertex curvePt2)
    {
        this.tangent = tangent;
        this.mass = mass;
        this.point = point;
        this.curvePt = curvePt;
        this.curvePt2 = curvePt2;
    }

    public TPEPointType GetTPEPointType()
    {
        if (curvePt != null)
            return TPEPointType.Point;
        else
            return TPEPointType.Cluster;
    }
}