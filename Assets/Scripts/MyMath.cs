using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class MyMath
{

    /// <summary>
    /// If intersecting, returns 0 <= t, u <= 1. Works in 2D with xz
    /// </summary>
    public static (float t, float u) Intersection(List<Vector3> points, int a_refPoint, int a_dirPoint, int b_refPoint, int b_dirPoint)
    {
        return Intersection(
            points[a_refPoint],
            points[a_dirPoint] - points[a_refPoint],
            points[b_refPoint],
            points[b_dirPoint] - points[b_refPoint]
        );
    }

    public static (float t, float u) Intersection(Vector3 a_refPoint, Vector3 a_dirWithLength, Vector3 b_refPoint, Vector3 b_dirWithLength)
    {
        var A = Matrix<float>.Build.DenseOfArray(new float[,]
        {
                    { a_dirWithLength.x, b_dirWithLength.x },
                    { a_dirWithLength.y, b_dirWithLength.y }
        });
        var b = Vector<float>.Build.Dense(new float[]
        {
                    b_refPoint.x - a_refPoint.x,
                    b_refPoint.y - a_refPoint.y
        });
        var values = A.Solve(b); // The values matrix is [t, -u]

        float t = values[0];
        float u = -values[1];

        return (t, u);
    }

    public static Vector3 xzToX0Z(Vector3 xz)
    {
        return new Vector3(xz.x, 0, xz.y);
    }

    public static Vector3 xyzToXZ(Vector3 xyz)
    {
        return new Vector3(xyz.x, xyz.z);
    }
}
