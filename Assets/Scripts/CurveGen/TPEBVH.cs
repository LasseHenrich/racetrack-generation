using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TPEBVH
{
    static TPEBVH instance;
    public static TPEBVH GetInstance { get { if (instance == null) instance = new TPEBVH(); return instance; } }

    public BVHNode2D CreateBVHFromCurve(EnergyCurve curve)
    {
        int numVerts = curve.NumVerts();
        List<VertexBody4D> verts = new List<VertexBody4D>(new VertexBody4D[numVerts]);

        // Loop over all the vertices
        for (int i = 0; i < numVerts; i++)
        {
            VertexBody4D currBody = VertToBody(curve, i);
            verts[currBody.elementIndex] = currBody;
        }

        BVHNode2D tree = new BVHNode2D(verts, 3, null, true);
        tree.RecomputeCentersOfMass(curve);
        BVHNode2D.globalID = 0;
        tree.RecursivelyAssignIDs();

        tree.numNodes = BVHNode2D.globalID;

        return tree;
    }

    private VertexBody4D VertToBody(EnergyCurve curve, int i)
    {
        CurveVertex p = curve.verts[i];
        Vector2 pos = p.Position();
        Vector2 tangent = p.Tangent();

        PosTan ptan = new PosTan(pos, tangent);
        float mass = p.AvgLength();
        int globalIndex = p.GlobalIndex();

        return new VertexBody4D(ptan, mass, globalIndex, BodyType.Vertex);
    }
}
