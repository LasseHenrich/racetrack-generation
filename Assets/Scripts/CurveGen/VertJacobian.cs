using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VertJacobian
{
    public Vector2 directional_x;
    public Vector2 directional_y;

    /// <summary>
    /// Interpret v as a row vector multiplied on the left.
    /// Then each entry is just the dot of v with the corresponding column.
    /// </summary>
    /// <param name="v"></param>
    /// <returns></returns>
    public Vector2 LeftMultiply(Vector2 v)
    {
        float x = Vector2.Dot(v, directional_x);
        float y = Vector2.Dot(v, directional_y);
        return new Vector2(x, y);
    }

    public static VertJacobian operator +(VertJacobian a, VertJacobian b)
    {
        Vector2 dir_x = a.directional_x + b.directional_x;
        Vector2 dir_y = a.directional_y + b.directional_y;
        return new VertJacobian { directional_x = dir_x, directional_y = dir_y };
    }

    public static VertJacobian operator -(VertJacobian a, VertJacobian b)
    {
        return a + (b * -1);
    }

    public static VertJacobian operator *(VertJacobian a, float c)
    {
        return new VertJacobian { directional_x = a.directional_x * c, directional_y = a.directional_y * c };
    }

    public static VertJacobian operator *(float c, VertJacobian a)
    {
        return a * c;
    }

    public static VertJacobian operator /(VertJacobian a, float c)
    {
        return a * (1 / c);
    }

    public static VertJacobian OuterProductToJacobian(Vector2 v1, Vector2 v2)
    {
        Vector2 col1 = v1 * v2.x;
        Vector2 col2 = v1 * v2.y;
        return new VertJacobian { directional_x = col1, directional_y = col2 };
    }

    internal void FillMatrix(Matrix<float> matrix)
    {
        float[] values = { directional_x.x, directional_y.x,
                           directional_x.y, directional_y.y};

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                matrix[i, j] = values[2 * i + j];
    }
}