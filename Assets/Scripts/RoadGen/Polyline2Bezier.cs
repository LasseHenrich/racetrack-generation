using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics.RootFinding;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Security.Cryptography;
using UnityEngine;

public static class PolylineToBezier
{
    static List<Vector3> bezierSpline;
    static float epsilon = 0.2f; // User-specified error that must be satisfied
    static float psi = 2f;       // In range [epsilon, psi) we try to minimize the error by reparameterization
    const int maxIterations = 10; // Reparameterization is tried maxIterations times 

    /// <summary>
    /// Converts a Polyline to a Bezier Spline
    /// </summary>
    /// <param name="d">Polyline</param>
    /// <param name="epsilon">User-specified error that must be satisfied</param>
    /// <param name="psi">In range [epsilon, psi) we try to minimize the error by reparameterization</param>
    /// <returns>A list of control- and anchor-points</returns>
    public static List<Vector3> Convert(List<Vector3> d, float epsilon, float psi)
    {
        PolylineToBezier.epsilon = epsilon;
        PolylineToBezier.psi = psi;

        bezierSpline = new List<Vector3>();
        d.Add(d[0]); // Add first point as last point

        Vector3 tHat1 = (d[1] - d[0]).normalized;   // Unit tangent 1
        Vector3 tHat2 = (d[^2] - d[^1]).normalized; // Unit tangent 2
        FitCubic(d, 0, d.Count - 1, tHat1, tHat2);

        bezierSpline.RemoveAt(bezierSpline.Count - 1); // Remove last anchor, as it's the same as first

        ResolveContinuity();

        return bezierSpline;
    }

    static void AddSegment(Vector3[] bezierCurve)
    {
        //bezierCurve.ToList().ForEach(x => Debug.Log(x));
        if (bezierSpline.Count == 0)
            bezierSpline.AddRange(bezierCurve);
        else
            bezierSpline.AddRange(bezierCurve.Skip(1));
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="d">Array of digitized points</param>
    /// <param name="first">Index of first point in region</param>
    /// <param name="last">Index of last point in region</param>
    /// <param name="tHat1">Endpoint 1 tangent</param>
    /// <param name="tHat2">Endpoint 2 tangent</param>
    static void FitCubic(List<Vector3> d, int first, int last, Vector3 tHat1, Vector3 tHat2)
    {
        //Debug.Log("FitCubic from " + first + " to " + last + " with tHat1 " + tHat1 + " and tHat2 " + tHat2);

        Vector3[] bezierCurve = new Vector3[4];
        int numPoints = last - first + 1;

        // If region only consists of two points, use heuristic
        if (numPoints == 2)
        {
            float dist = Vector3.Distance(d[last], d[first]) / 3f;

            bezierCurve[0] = d[first];
            bezierCurve[3] = d[last];
            bezierCurve[1] = bezierCurve[0] + tHat1 * dist;
            bezierCurve[2] = bezierCurve[3] + tHat2 * dist;
            AddSegment(bezierCurve);
            return;
        }


        // Parameterize points and attempt to fit curve
        List<float> u = ChordLengthParameterization(d, first, last);
        //u.ForEach(x => Debug.Log(x));
        bezierCurve = GenerateBezier(d, first, last, u, tHat1, tHat2);

        //bezierCurve.ToList().ForEach(x => Debug.Log(x));

        // Compute the maximum distance from points to curve
        int splitPoint = -1;
        float error = ComputeMaxError(d, first, last, bezierCurve, u, ref splitPoint);
        //Debug.Log(error);

        if (error < epsilon)
        {
            //Debug.Log("Error less than epsilon " + epsilon);
            AddSegment(bezierCurve);
            return;
        }

        // If the error is not too large, try some reparameterization
        if (error < psi)
        {
            //Debug.Log("Error less than psi");
            for (int i = 0; i < maxIterations; i++)
            {
                List<float> uPrime = Reparameterize(d, first, last, u, bezierCurve);
                bezierCurve = GenerateBezier(d, first, last, uPrime, tHat1, tHat2);
                error = ComputeMaxError(d, first, last, bezierCurve, uPrime, ref splitPoint);
                if (error < epsilon)
                {
                    AddSegment(bezierCurve);
                    return;
                }
                u = uPrime;
            }
            Debug.Log("maxIterations " + maxIterations + " exceeded");
        }

        // Fitting failed -- Split at maximum error point and fit recursively
        Vector3 tHatCenter = ComputeCenterTangent(d, splitPoint);
        FitCubic(d, first, splitPoint, tHat1, tHatCenter);
        FitCubic(d, splitPoint, last, -tHatCenter, tHat2);
    }

    /// <summary>
    /// REturns the tangent of "center" that is part of "d"
    /// </summary>
    /// <param name="d">Digitized points</param>
    /// <param name="center">Index of point inside region</param>
    /// <returns></returns>
    static Vector3 ComputeCenterTangent(List<Vector3> d, int center)
    {
        Vector3 V1 = d[center - 1] - d[center];
        Vector3 V2 = d[center] - d[center + 1];
        return (V1 + V2).normalized;
    }

    /// <summary>
    /// Given a set of points d and their parameterization u, this function tries to find a better parameterization
    /// </summary>
    /// <param name="d">Aary of digitized points</param>
    /// <param name="first">Index of first point in region</param>
    /// <param name="last">Index of last point in region</param>
    /// <param name="u">Current parameterization</param>
    /// <param name="bezierCurve">Current fitted curve</param>
    /// <returns></returns>
    static List<float> Reparameterize(List<Vector3> d, int first, int last, List<float> u, Vector3[] bezierCurve)
    {
        List<float> uPrime = new();
        Vector3[] _bezierCurve = (Vector3[])bezierCurve.Clone();

        for (int i = first; i <= last; i++)
            uPrime.Add(NewtonRaphsonRootFind(_bezierCurve, d[i], u[i - first]));
        return uPrime;
    }

    /// <summary>
    /// Uses the Newon-Raphson iteration to find better root
    /// </summary>
    /// <param name="Q">Current fitted curve</param>
    /// <param name="P">Digitized point</param>
    /// <param name="u">Parameter value for "P"</param>
    /// <returns></returns>
    static float NewtonRaphsonRootFind(Vector3[] Q, Vector3 P, float u)
    {
        //List<Vector2> Q = new(_Q);

        // Compute Q(u)
        Vector3 Q_u = Bezier(Q, u);

        // Generate control vetices for Q'
        Vector3[] Q1 = new Vector3[3];
        for (int i = 0; i < 3; i++)
            Q1[i] = (Q[i + 1] - Q[i]) * 3f;

        // Geneate control vertices for Q''
        Vector3[] Q2 = new Vector3[2];
        for (int i = 0; i < 2; i++)
            Q2[i] = (Q1[i + 1] - Q1[i]) * 2f;

        // Compute Q'(u) and Q''(u)
        Vector3 Q1_u = Bezier(Q1, u);
        Vector3 Q2_u = Bezier(Q2, u);

        // Comptue f(u)/f'(u)
        float numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
        float denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
                            (Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);
        return denominator == 0f ? u : (u - (numerator / denominator));
    }

    /// <summary>
    /// Evaluates a Bezier curve defined by V with degree V.Count at a particulat parameter value
    /// </summary>
    /// <param name="V"></param>
    /// <param name="t"></param>
    static Vector3 Bezier(Vector3[] V, float t)
    {
        int degree = V.Length - 1;
        List<Vector3> Vtemp = new(V);

        // Triangle computation
        for (int i = 1; i <= degree; i++)
            for (int j = 0; j <= degree - i; j++)
                Vtemp[j] = (1f - t) * Vtemp[j] + t * Vtemp[j + 1];

        return Vtemp[0];
    }

    /// <summary>
    /// Finds the maximum squared distance of the digitized points d to the fitted Curve bezierCurve
    /// </summary>
    /// <param name="d">Array of digitized points</param>
    /// <param name="first">Index of first point in region</param>
    /// <param name="last">Index of last poitn in region</param>
    /// <param name="bezierCurve">Fitted Bezier curve</param>
    /// <param name="u">Parameterization of points</param>
    /// <param name="splitPoint">Point of maximum error</param>
    /// <returns></returns>
    static float ComputeMaxError(List<Vector3> d, int first, int last, Vector3[] bezierCurve, List<float> u, ref int splitPoint)
    {
        splitPoint = (last - first + 1) / 2;
        float maxDist = 0f;
        for (int i = first + 1; i < last; i++)
        {
            Vector3 P = Bezier(bezierCurve, u[i - first]);
            Vector3 V = P - d[i];
            float dist = V.sqrMagnitude;
            if (dist >= maxDist)
            {
                maxDist = dist;
                splitPoint = i;
            }
        }

        return maxDist;
    }

    /// <summary>
    /// Uses least-squares method to find Bezier control points for region
    /// </summary>
    /// <param name="d">Array of digitized points</param>
    /// <param name="first">Index of first point in region</param>
    /// <param name="last">Index of last point in region</param>
    /// <param name="uPrime">Parameterization of region</param>
    /// <param name="tHat1">Unit tangent at endpoint 1</param>
    /// <param name="tHat2">Unit tangent at endpoint 2</param>
    /// <returns></returns>
    static Vector3[] GenerateBezier(List<Vector3> d, int first, int last, List<float> uPrime, Vector3 tHat1, Vector3 tHat2)
    {
        Vector3[] bezierCurve = new Vector3[4];
        int numPoints = last - first + 1;

        // Compute A
        Vector3[,] A = new Vector3[numPoints, 2]; // Percomputed rhs of equation
        for (int i = 0; i < numPoints; i++)
        {
            A[i, 0] = tHat1 * B1(uPrime[i]);
            A[i, 1] = tHat2 * B2(uPrime[i]);
        }

        // Create the C and X matrices
        float[,] C = new float[2, 2];
        float[] X = new float[2];

        for (int i = 0; i < numPoints; i++)
        {
            C[0, 0] += Vector3.Dot(A[i, 0], A[i, 0]);
            C[0, 1] += Vector3.Dot(A[i, 0], A[i, 1]);
            C[1, 0] = C[0, 1];
            C[1, 1] += Vector3.Dot(A[i, 1], A[i, 1]);
            //Debug.Log(C[0, 0]);
            //Debug.Log(C[0, 1]);
            //Debug.Log(C[1, 0]);
            //Debug.Log(C[1, 1]);

            Vector3 tmp = d[first + i] - (
                d[first] * B0(uPrime[i]) +
                d[first] * B1(uPrime[i]) +
                d[last] * B2(uPrime[i]) +
                d[last] * B3(uPrime[i])
                );

            X[0] += Vector3.Dot(A[i, 0], tmp);
            X[1] += Vector3.Dot(A[i, 1], tmp);
        }

        // Compute the determinants of C and X
        float det_C0_C1 = C[0, 0] * C[1, 1] - C[1, 0] * C[0, 1];
        float det_C0_X = C[0, 0] * X[1] - C[1, 0] * X[0];
        float det_X_C1 = X[0] * C[1, 1] - X[1] * C[0, 1];
        //Debug.Log(det_C0_C1);
        //Debug.Log(det_C0_X);
        //Debug.Log(det_X_C1);

        // Derive alpha values
        float alpha_l = det_C0_C1 == 0 ? 0f : det_X_C1 / det_C0_C1;
        float alpha_r = det_C0_C1 == 0 ? 0f : det_C0_X / det_C0_C1;
        //Debug.Log(alpha_l);
        //Debug.Log(alpha_r);

        // If alpha is negative, use the Wu/Barsky heuristic.
        // alpha = 0 would lead to coincident control points, resulting in a
        // division by zero in any subsequent NewtonRaphsonRootFind()-call.
        float segLength = Vector2.Distance(d[last], d[first]);
        float epsilon = 1e-6f * segLength;
        if (alpha_l < epsilon || alpha_r < epsilon)
        {
            //Debug.Log("alpha close to zero");
            // Fall back on standard (probably inaccurate) formula and subdivide further if needed
            float dist = segLength / 3f;
            bezierCurve[0] = d[first];
            bezierCurve[3] = d[last];
            bezierCurve[1] = bezierCurve[0] + tHat1 * dist;
            bezierCurve[2] = bezierCurve[3] + tHat2 * dist;
            return bezierCurve;
        }

        // First and last control points of the Bezier curve are
        // positioned exactly at the fist and last data points.
        // Control points 1 and 2 are positioned an alpha distance out
        // on the tangent vectors, left and right, respectively
        bezierCurve[0] = d[first];
        bezierCurve[3] = d[last];
        bezierCurve[1] = bezierCurve[0] + tHat1 * alpha_l;
        bezierCurve[2] = bezierCurve[3] + tHat2 * alpha_r;
        return bezierCurve;
    }

    static void ResolveContinuity()
    {
        Vector3 firstControlVec = bezierSpline[1] - bezierSpline[0];
        float firstControlLength = firstControlVec.magnitude;
        Vector3 lastControlVec = bezierSpline[^1] - bezierSpline[0];
        float lastControlLength = lastControlVec.magnitude;

        Vector3 newFirstControlVec = (firstControlVec * firstControlLength - lastControlVec * lastControlLength).normalized * firstControlLength; // Weighted Average
        Vector3 newLastControlVec = (lastControlVec * lastControlLength - firstControlVec * firstControlLength).normalized * lastControlLength;

        bezierSpline[1] = bezierSpline[0] + newFirstControlVec;
        bezierSpline[^1] = bezierSpline[0] + newLastControlVec;
    }

    static float B0(float u)
    {
        float tmp = 1f - u;
        return tmp * tmp * tmp;
    }

    static float B1(float u)
    {
        float tmp = 1f - u;
        return 3 * u * tmp * tmp;
    }

    static float B2(float u)
    {
        float tmp = 1f - u;
        return 3 * u * u * tmp;
    }
    static float B3(float u)
    {
        return u * u * u;
    }

    /// <summary>
    /// Assign parameter values to digized points using relative distances between points
    /// </summary>
    /// <param name="d">Array of digitized points</param>
    /// <param name="first">Index of first point in region</param>
    /// <param name="last">Index of last point in region</param>
    /// <returns></returns>
    static List<float> ChordLengthParameterization(List<Vector3> d, int first, int last)
    {
        List<float> u = new()
        {
            0
        };

        for (int i = first + 1; i <= last; i++)
            u.Add(u[i - first - 1] + Vector3.Distance(d[i], d[i - 1]));

        for (int i = first + 1; i <= last; i++)
            u[i - first] /= u[last - first];

        return u;
    }
}
