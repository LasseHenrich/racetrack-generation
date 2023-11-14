using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CurveGenUtils
{
    //static CurveGenUtils instance;
    //public static CurveGenUtils Instance { get { if (instance == null) instance = new CurveGenUtils(); return instance; } }

    public static void AddToRow(Matrix<float> A, int row, Vector2 toAdd)
    {
        A[row, 0] += toAdd.x;
        A[row, 1] += toAdd.y;
    }

    public static Vector2 SelectRow(Matrix<float> A, int row)
    {
        return new Vector2(A[row, 0], A[row, 1]);
    }

    public static void SetRow(Matrix<float> A, int row, Vector2 toAdd)
    {
        A[row, 0] = toAdd.x;
        A[row, 1] = toAdd.y;
    }

    internal static Vector2 VectorMax(Vector2 a, Vector2 b)
    {
        return new Vector2(Mathf.Max(a.x, b.x), Mathf.Max(a.y, b.y));
    }

    internal static Vector2 VectorMin(Vector2 a, Vector2 b)
    {
        return new Vector2(Mathf.Min(a.x, b.x), Mathf.Min(a.y, b.y));
    }

    internal static void CrossMatrix(Matrix<float> skw)
    {
        skw[0, 0] = 0; skw[0, 1] = 1;
        skw[1, 0] = -1; skw[1, 1] = 0;
    }
    internal static Matrix<float> CrossMatrix()
    {
        Matrix<float> skw = Matrix<float>.Build.Dense(2, 2);
        skw[0, 0] = 0; skw[0, 1] = 1;
        skw[1, 0] = -1; skw[1, 1] = 0;
        return skw;
    }

    internal static void CrossMatrix3D(Vector3 v, Matrix<float> skw)
    {
        skw[0, 0] = 0; skw[0, 1] = -v.z; skw[0, 2] = v.y;
        skw[1, 0] = v.z; skw[1, 1] = 0; skw[1, 2] = -v.x;
        skw[2, 0] = -v.y; skw[2, 1] = v.x; skw[2, 2] = 0;
    }


    internal static float Cross(Vector2 v1, Vector2 v2)
    {
        return (v1.x * v2.y) - (v1.y * v2.x);
    }
}
