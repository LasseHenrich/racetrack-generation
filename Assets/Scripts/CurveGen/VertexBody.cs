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
    public Vector2 position;
    public Vector2 tangent;
    public float mass;
    public int vertIndex;
}

public class PosTan
{
    public Vector2 position;
    public Vector2 tangent;

    public PosTan(Vector2 position, Vector2 tangent)
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

public class VertexBody4D
{
    public PosTan pt;
    public float mass;
    public int elementIndex;
    public BodyType type;

    public VertexBody4D(PosTan pt, float mass, int elementIndex, BodyType type)
    {
        this.pt = pt;
        this.mass = mass;
        this.elementIndex = elementIndex;
        this.type = type;
    }

    public VertexBody4D()
    {
        this.pt = new PosTan();
    }
}