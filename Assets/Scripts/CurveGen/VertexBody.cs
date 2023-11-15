using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum BodyType
{
    Vertex,
    Edge
}

public class VertexBody
{
    public Vector3 position;
    public Vector3 tangent;
    public float mass;
    public int vertIndex;
}

public class PosTan
{
    public Vector3 position;
    public Vector3 tangent;

    public PosTan(Vector3 position, Vector3 tangent)
    {
        this.position = position;
        this.tangent = tangent;
    }

    public PosTan()
    {

    }

    public void Print()
    {
        Debug.Log(position + "; " + tangent);
    }
}

public class VertexBody6D
{
    public PosTan pt;
    public float mass;
    public int elementIndex;
    public BodyType type;

    public VertexBody6D(PosTan pt, float mass, int elementIndex, BodyType type)
    {
        this.pt = pt;
        this.mass = mass;
        this.elementIndex = elementIndex;
        this.type = type;
    }

    public VertexBody6D()
    {
        this.pt = new PosTan();
    }
}