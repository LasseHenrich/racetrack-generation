using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using Unity.Rendering.HybridV2;
using UnityEngine;

public abstract class Potential
{
    protected readonly EnergyCurve curve;
    protected readonly float weight;
    /// <summary>
    /// If true, this potential will only be active when the target length is reached
    /// </summary>
    public readonly bool delayed;

    public Potential(EnergyCurve curve, float weight, bool delayed)
    {
        this.curve = curve;
        this.weight = weight;
        this.delayed = delayed;
    }

    public abstract float CurrentValue();
    public abstract void AddGradient(Matrix<float> gradient);
}

public class VectorFieldPotential : Potential
{

    public VectorFieldPotential(EnergyCurve curve, float weight, bool delayed) : base(curve, weight, delayed)
    {
    }

    // Todo: Have a look at this again -- seems different than https://github.com/icethrush/repulsive-curves/blob/master/src/extra_potentials.cpp
    // and should maybe work differently in 3D, even for our simpler "VectorField"
    public override void AddGradient(Matrix<float> gradient)
    {
        int numEdges = curve.NumEdges();

        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge e_i = curve.edges[i];

            Vector3 tangent = e_i.Tangent();

            // Sample the vector field at edge midpoint
            Vector3 vf = Sample(tangent);
            Vector3 vxt = Vector3.Cross(vf, tangent);

            Matrix<float> Y_skw = Matrix<float>.Build.Dense(3, 3);
            CurveGenUtils.CrossMatrix3D(vf, Y_skw);

            CurveVertex v_prev = e_i.GetPrevVertex();
            CurveVertex v_next = e_i.GetNextVertex();

            // Derivatives of lengths and tangents
            VertJacobian ddx_t_prev_J = TPE.GetInstance.EdgeTangentWrtVert(e_i, v_prev); // Derivative of Tangent wrt Prev, as Jacobian
            VertJacobian ddx_t_next_J = TPE.GetInstance.EdgeTangentWrtVert(e_i, v_next);
            Vector3 deriv_len_prev = TPE.GetInstance.EdgeLengthWrtVert(e_i, v_prev); // Derivative of Length wrt Prev
            Vector3 deriv_len_next = TPE.GetInstance.EdgeLengthWrtVert(e_i, v_next);

            // Transfer tangent Jacobians into 3x3 matrices
            Matrix<float> ddx_t_prev = Matrix<float>.Build.Dense(3, 3), ddx_t_next = Matrix<float>.Build.Dense(3, 3);
            ddx_t_prev_J.FillMatrix(ddx_t_prev); // Derivative of Tangent wrt Prev
            ddx_t_next_J.FillMatrix(ddx_t_next);

            // Transfer Unity-Vector into MathNet-Vector
            Vector<float> vf_v = Vector<float>.Build.Dense(3);
            vf_v[0] = vf.x; vf_v[1] = vf.y; vf_v[2] = vf.z;
            Vector<float> vxt_v = Vector<float>.Build.Dense(3);
            vxt_v[0] = vxt.x; vxt_v[1] = vxt.y; vxt_v[2] = vxt.z;

            Vector<float> deriv_prev = -2 * (Y_skw * ddx_t_prev) * vxt_v;
            Vector<float> deriv_next = -2 * (Y_skw * ddx_t_next) * vxt_v;
            float normcross2 = vxt_v * vxt_v;
            float len = e_i.Length();

            Vector3 add_prev = ((deriv_len_prev * normcross2) + (len * new Vector3(deriv_prev[0], deriv_prev[1], deriv_prev[3]))) * weight;
            Vector3 add_next = ((deriv_len_next * normcross2) + (len * new Vector3(deriv_next[0], deriv_next[1], deriv_next[3]))) * weight;
            CurveGenUtils.AddToRow(gradient, v_prev.GlobalIndex(), add_prev);
            CurveGenUtils.AddToRow(gradient, v_next.GlobalIndex(), add_next);
        }
    }

    public override float CurrentValue()
    {
        int numEdges = curve.NumEdges();
        float sum = 0;

        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge e_i = curve.edges[i];
            Vector3 tangent = e_i.Tangent();
            Vector3 vf = Sample(tangent);
            Vector3 vxt = Vector3.Cross(vf, tangent);
            sum += Vector3.Dot(vxt, vxt);
        }

        return sum * weight;
    }

    public Vector3 Sample(Vector3 tangent)
    {
        //float angle = Mathf.Atan2(Mathf.Abs(tangent.y), Mathf.Abs(tangent.x));
        float angle = Mathf.Atan(tangent.y / tangent.x);
        //Debug.Log(angle);
        return angle >= 0 ? new Vector3(1, 1, 0) : new Vector3(-1, 1, 0);
    }

    public VertJacobian SpatialDerivative()
    {
        return new VertJacobian { directional_x = Vector3.zero, directional_y = Vector3.zero };
    }
}

public enum PotentialType
{
    VectorField
}