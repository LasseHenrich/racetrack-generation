using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SoboSlobo
{
    static SoboSlobo instance;
    public static SoboSlobo GetInstance { get { if (instance == null) instance = new SoboSlobo(); return instance; } }

    float alpha { get { return EnergyCurve.alpha; } }
    float beta { get { return EnergyCurve.beta; } }

    public void Sobolev3XWithConstraints(EnergyCurve curve, ref Matrix<float> sps1, float diagEps = 0)
    {
        int numRows = curve.constraintSet.NumConstraintRows() + curve.constraintSet.NumExpectedCols();
        sps1 = Matrix<float>.Build.Dense(numRows, numRows);
        //A = Matrix<float>.Build.Dense(network.NumVerts() * 2 + 2, network.NumVerts() * 2 + 2); // DEBUG
        SobolevGramMatrix3X(curve, sps1, diagEps);
        curve.constraintSet.FillDenseBlock(sps1);
    }

    void SobolevGramMatrix3X(EnergyCurve curve, Matrix<float> sps1, float diagEps)
    {
        int numVerts = curve.NumVerts();
        Matrix<float> topLeft = Matrix<float>.Build.Dense(numVerts, numVerts); // This is A_line

        SobolevGramMatrix(curve, topLeft, diagEps);

        for (int i = 0; i < numVerts; i++)
        {
            for (int j = 0; j < numVerts; j++)
            {
                sps1[3 * i + 0, 3 * j + 0] = topLeft[i, j];
                sps1[3 * i + 1, 3 * j + 1] = topLeft[i, j];
                sps1[3 * i + 2, 3 * j + 2] = topLeft[i, j];
            }
        }
    }

    void SobolevGramMatrix(EnergyCurve curve, Matrix<float> A, float diagEps)
    {
        int numVerts = curve.NumVerts();
        int numEdges = curve.NumEdges();

        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge pc_i = curve.edges[i];

            for (int j = 0; j < numEdges; j++)
            {
                CurveEdge pc_j = curve.edges[j];
                if (pc_i == pc_j || pc_i.IsNeighbors(pc_j)) continue;

                AddEdgePairContribution(pc_i, pc_j, A);
                AddEdgePairContributionLow(pc_i, pc_j, A);
            }
        }

        for (int i = 0; i < numVerts; i++)
        {
            A[i, i] += diagEps * curve.verts[i].AvgLength();
        }
    }

    void AddEdgePairContribution(CurveEdge s, CurveEdge t, Matrix<float> A)
    {
        List<CurveVertex> endpoints = new() { s.GetPrevVertex(), s.GetNextVertex(), t.GetPrevVertex(), t.GetNextVertex() };
        float len_s = s.Length();
        float len_r = t.Length();
        Vector3 mid_s = s.Midpoint();
        Vector3 mid_t = t.Midpoint();
        Vector3 tangent_s = s.Tangent();
        Vector3 tangent_t = t.Tangent();

        float dist_term = MetricDistanceTerm(mid_s, mid_t, tangent_s, tangent_t);

        foreach (CurveVertex u in endpoints)
        {
            foreach (CurveVertex v in endpoints)
            {
                Vector3 u_hat_s = HatGradientOnEdge(s, u);
                Vector3 u_hat_t = HatGradientOnEdge(t, u);
                Vector3 v_hat_s = HatGradientOnEdge(s, v);
                Vector3 v_hat_t = HatGradientOnEdge(t, v);

                float numer = Vector3.Dot(u_hat_s - u_hat_t, v_hat_s - v_hat_t);
                int index_u = u.GlobalIndex();
                int index_v = v.GlobalIndex();

                A[index_u, index_v] += numer * dist_term * len_s * len_r;
            }
        }
    }

    void AddEdgePairContributionLow(CurveEdge s, CurveEdge t, Matrix<float> A)
    {
        List<CurveVertex> endpoints = new() { s.GetPrevVertex(), s.GetNextVertex(), t.GetPrevVertex(), t.GetNextVertex() };
        float len_1 = s.Length();
        float len_2 = t.Length();
        Vector3 mid_s = s.Midpoint();
        Vector3 mid_t = t.Midpoint();
        Vector3 tangent_s = s.Tangent();
        Vector3 tangent_t = t.Tangent();

        float kf_st = MetricDistanceTermLow(mid_s, mid_t, tangent_s, tangent_t);

        foreach (CurveVertex u in endpoints)
        {
            foreach (CurveVertex v in endpoints)
            {
                float u_s = HatMidpoint(s, u);
                float u_t = HatMidpoint(t, u);
                float v_s = HatMidpoint(s, v);
                float v_t = HatMidpoint(t, v);

                float numer = (u_s - u_t) * (v_s - v_t);
                int index_u = u.GlobalIndex();
                int index_v = v.GlobalIndex();

                A[index_u, index_v] += numer * kf_st * len_1 * len_2;
            }
        }
    }

    float MetricDistanceTerm(Vector3 v1, Vector3 v2, Vector3 t1, Vector3 t2)
    {
        float s_pow = (beta - 1) / alpha;
        float dist_term = 1.0f / Mathf.Pow((v1 - v2).magnitude, 2 * (s_pow - 1) + 1);
        return dist_term;
    }

    float MetricDistanceTermLow(Vector3 v1, Vector3 v2, Vector3 t1, Vector3 t2)
    {
        float s_pow = (beta - 1) / alpha;
        s_pow = 2 * (s_pow - 1) + 1;
        float a = 2f;
        float b = 4f + s_pow;
        return TPE.GetInstance.TPE_Kf_pts_sym(v1, v2, t1, t2, a, b);
    }

    Vector3 HatGradientOnEdge(CurveEdge edge, CurveVertex vertex)
    {
        CurveVertex edgeStart = edge.GetPrevVertex();
        CurveVertex edgeEnd = edge.GetNextVertex();

        if (vertex.Equals(edgeStart))
        {
            Vector3 towardsVertex = edgeStart.Position() - edgeEnd.Position();
            float length = towardsVertex.magnitude;
            // 1 over length times normalized edge vector towards the vertex
            return towardsVertex / (length * length);
        }
        else if (vertex.Equals(edgeEnd))
        {
            Vector3 towardsVertex = edgeEnd.Position() - edgeStart.Position();
            float length = towardsVertex.magnitude;
            return towardsVertex / (length * length);
        }
        else
        {
            return Vector3.zero;
        }
    }

    float HatMidpoint(CurveEdge edge, CurveVertex vertex)
    {
        CurveVertex edgeStart = edge.GetPrevVertex();
        CurveVertex edgeEnd = edge.GetNextVertex();

        if (vertex.Equals(edgeStart))
            return 0.5f;
        else if (vertex.Equals(edgeEnd))
            return 0.5f;
        else
            return 0f;
    }
}
