using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TPEBVH
{
    static TPEBVH instance;
    public static TPEBVH GetInstance { get { if (instance == null) instance = new TPEBVH(); return instance; } }

    public BVHNode3D CreateBVHFromCurve(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        List<VertexBody6D> verts = new List<VertexBody6D>(new VertexBody6D[numVerts]);

        // Loop over all the vertices
        for (int i = 0; i < numVerts; i++)
        {
            VertexBody6D currBody = VertToBody(curve, i);
            verts[currBody.elementIndex] = currBody;
        }

        BVHNode3D tree = new BVHNode3D(verts, 3, null, true);
        tree.RecomputeCentersOfMass(curve);
        BVHNode3D.globalID = 0;
        tree.RecursivelyAssignIDs();

        tree.numNodes = BVHNode3D.globalID;

        return tree;
    }

    private VertexBody6D VertToBody(EnergyCurve curve, int i)
    {
        CurveVertex p = curve.verts[i];
        Vector2 pos = p.Position();
        Vector2 tangent = p.Tangent();

        PosTan ptan = new PosTan(pos, tangent);
        float mass = p.AvgLength();
        int globalIndex = p.GlobalIndex();

        return new VertexBody6D(ptan, mass, globalIndex, BodyType.Vertex);
    }
}
