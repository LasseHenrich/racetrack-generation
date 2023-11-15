using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CurveGenUtils
{
    //static CurveGenUtils instance;
    //public static CurveGenUtils Instance { get { if (instance == null) instance = new CurveGenUtils(); return instance; } }

    public static void AddToRow(Matrix<float> A, int row, Vector3 toAdd)
    {
        A[row, 0] += toAdd.x;
        A[row, 1] += toAdd.y;
        A[row, 2] += toAdd.z;
    }

    public static Vector3 SelectRow(Matrix<float> A, int row)
    {
        return new Vector3(A[row, 0], A[row, 1], A[row, 2]);
    }

    public static void SetRow(Matrix<float> A, int row, Vector3 toAdd)
    {
        A[row, 0] = toAdd.x;
        A[row, 1] = toAdd.y;
        A[row, 2] = toAdd.z;
    }

    internal static Vector3 VectorMax(Vector3 a, Vector3 b)
    {
        return new Vector3(Mathf.Max(a.x, b.x), Mathf.Max(a.y, b.y), Mathf.Max(a.z, b.z));
    }

    internal static Vector3 VectorMin(Vector3 a, Vector3 b)
    {
        return new Vector3(Mathf.Min(a.x, b.x), Mathf.Min(a.y, b.y), Mathf.Min(a.z, b.z));
    }

    /*
    internal static void CrossMatrix2D(Matrix<float> skw)
    {
        skw[0, 0] = 0; skw[0, 1] = 1;
        skw[1, 0] = -1; skw[1, 1] = 0;
    }
    */

    /*
    internal static Matrix<float> CrossMatrix2D()
    {
        Matrix<float> skw = Matrix<float>.Build.Dense(2, 2);
        skw[0, 0] = 0; skw[0, 1] = 1;
        skw[1, 0] = -1; skw[1, 1] = 0;
        return skw;
    }
    */

    internal static void CrossMatrix3D(Vector3 v, Matrix<float> skw)
    {
        skw[0, 0] = 0; skw[0, 1] = -v.z; skw[0, 2] = v.y;
        skw[1, 0] = v.z; skw[1, 1] = 0; skw[1, 2] = -v.x;
        skw[2, 0] = -v.y; skw[2, 1] = v.x; skw[2, 2] = 0;
    }

    /*
    internal static float Cross2D(Vector2 v1, Vector2 v2)
    {
        return (v1.x * v2.y) - (v1.y * v2.x);
    }
    */
}
