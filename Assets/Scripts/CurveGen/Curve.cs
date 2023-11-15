using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

[Serializable]
public class Curve
{
    public List<CurveVertex> verts;
    public List<CurveEdge> edges;
    public Matrix<float> positions; // This is not serializable

    protected bool curveClosed;

    #region Serialization
    public List<Vector3> s_posList;
    #endregion

    /// <summary>
    /// Initializes Vertices and Edges from a position list. Note that the position list
    /// will be treated as a list of continuously connected vertices.
    /// </summary>
    /// <param name="posList"></param>
    protected void InitVertsEdgesFromPositions(List<Vector3> posList)
    {
        int numVerts = posList.Count;
        positions = Matrix<float>.Build.Dense(numVerts, 3);

        verts = new List<CurveVertex>();
        for (int i = 0; i < numVerts; i++)
        {
            verts.Add(new CurveVertex());
        }

        int numEdges = numVerts - 1;
        if (curveClosed) numEdges++;

        edges = new List<CurveEdge>();
        for (int i = 0; i < numEdges; i++) // For every but the last vertex...
        {
            edges.Add(new CurveEdge()); // Create a new edge to the next vertex
            edges[i].Init(verts[i], verts[(i + 1) % numVerts], i);
        }

        for (int i = 0; i < numVerts; i++)
        {
            verts[i].Init(this, edges[(i + numEdges - 1) % numEdges], edges[i], i);
            verts[i].SetPosition(posList[i]);
        }
    }

    public int NumVerts() => verts.Count;

    public int NumEdges() => edges.Count;

    public float TotalLength()
    {
        float length = 0;
        for (int i = 0; i < NumEdges(); i++)
            length += edges[i].Length();
        return length;
    }

    public void TryInitVertsEdgesFromSerialized()
    {
        InitVertsEdgesFromPositions(s_posList);
    }

    public void SerializeVertsEdgesPositions()
    {
        s_posList = new List<Vector3>();
        for (int i = 0; i < positions.RowCount; i++)
        {
            s_posList.Add(new(positions.Row(i)[0], positions.Row(i)[1]));
        }
    }
}

public class CurveVertex
{
    Curve network;
    CurveEdge prevEdge;
    CurveEdge nextEdge;
    int index;

    public void Init(Curve network, CurveEdge prevEdge, CurveEdge nextEdge, int index)
    {
        this.network = network;
        this.prevEdge = prevEdge;
        this.nextEdge = nextEdge;
        this.index = index;
    }

    public Vector3 Tangent()
    {
        Vector3 tangent = Vector3.zero;
        for (int i = 0; i < NumEdges(); i++)
            tangent += Edge(i).Tangent();
        return tangent.normalized;
    }

    public int GlobalIndex() => index;

    public int NumEdges() => 2;

    public CurveEdge Edge(int i) => i == 0 ? prevEdge : nextEdge;

    /// <summary>
    /// (DualLength) Returns the average length of all attacked edges
    /// </summary>
    /// <returns></returns>
    public float AvgLength()
    {
        float length = 0;
        for (int i = 0; i < NumEdges(); i++)
            length += Edge(i).Length();
        return length / NumEdges();
    }

    public Vector3 Position() => CurveGenUtils.SelectRow(network.positions, index);

    public void SetPosition(Vector3 newPos) => CurveGenUtils.SetRow(network.positions, index, newPos);
}

public class CurveEdge
{
    CurveVertex prevVertex;
    CurveVertex nextVertex;
    int index;

    public void Init(CurveVertex prevVertex, CurveVertex nextVertex, int index)
    {
        this.prevVertex = prevVertex;
        this.nextVertex = nextVertex;
        this.index = index;
    }

    public CurveVertex GetNextVertex() => nextVertex;
    public CurveVertex GetPrevVertex() => prevVertex;
    public int GetIndex() => index;

    public Vector3 Tangent() => (nextVertex.Position() - prevVertex.Position()).normalized;

    public float Length() => (nextVertex.Position() - prevVertex.Position()).magnitude;

    public CurveVertex Opposite(CurveVertex v)
    {
        if (nextVertex.Equals(v)) return prevVertex;
        else if (prevVertex.Equals(v)) return nextVertex;
        else throw new System.ArgumentException("v is wrong");
    }

    public bool IsNeighbors(CurveEdge other)
    {
        bool shared = (nextVertex == other.nextVertex) ||
            (nextVertex == other.prevVertex) ||
            (prevVertex == other.nextVertex) ||
            (prevVertex == other.prevVertex);
        return shared && (!other.Equals(this));
    }

    public Vector3 Midpoint() => (nextVertex.Position() + prevVertex.Position()) * 0.5f;
}