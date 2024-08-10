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

    // This function calculates the gradient for an EDGE,
    // and therefore TWO Vertices. However, the calculations
    // could be handled VERTEX-INDEPENDENT, with a
    // TWO-DIMENSIONAL Gradient instead.
    public override void AddGradient(Matrix<float> gradient)
    {
        int numEdges = curve.NumEdges();

        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge e_i = curve.edges[i];

            Vector2 tangent = e_i.Tangent();

            // Sample the vector field at edge midpoint
            Vector2 vf = Sample(tangent);
            float vxt = CurveGenUtils.Cross(vf, tangent);

            CurveVertex v_prev = e_i.GetPrevVertex();
            CurveVertex v_next = e_i.GetNextVertex();

            // Derivatives of lengths and tangents
            VertJacobian ddx_t_prev_J = TPE.GetInstance.EdgeTangentWrtVert(e_i, v_prev); // Derivative of Tangent wrt Prev, as Jacobian
            VertJacobian ddx_t_next_J = TPE.GetInstance.EdgeTangentWrtVert(e_i, v_next);
            Vector2 deriv_len_prev = TPE.GetInstance.EdgeLengthWrtVert(e_i, v_prev); // Derivative of Length wrt Prev
            Vector2 deriv_len_next = TPE.GetInstance.EdgeLengthWrtVert(e_i, v_next);

            // Transfer tangent Jacobians into 2x2 matrices
            Matrix<float> ddx_t_prev = Matrix<float>.Build.Dense(2, 2), ddx_t_next = Matrix<float>.Build.Dense(2, 2);
            ddx_t_prev_J.FillMatrix(ddx_t_prev); // Derivative of Tangent wrt Prev
            ddx_t_next_J.FillMatrix(ddx_t_next);

            // Transfer Unity-Vector into MathNet-Vector
            Vector<float> vf_v = Vector<float>.Build.Dense(2);
            vf_v[0] = vf.x; vf_v[1] = vf.y;

            Vector<float> deriv_prev = -2 * (vf_v * ddx_t_prev) * vxt;
            Vector<float> deriv_next = -2 * (vf_v * ddx_t_next) * vxt;
            float normcross2 = vxt * vxt;
            float len = e_i.Length();

            Vector2 add_prev = ((deriv_len_prev * normcross2) + (len * new Vector2(deriv_prev[0], deriv_prev[1]))) * weight;
            Vector2 add_next = ((deriv_len_next * normcross2) + (len * new Vector2(deriv_next[0], deriv_next[1]))) * weight;
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
            Vector2 tangent = e_i.Tangent();
            Vector2 vf = Sample(tangent);
            float vxt = CurveGenUtils.Cross(vf, tangent);
            sum += vxt * vxt;
        }

        return sum * weight;
    }

    public Vector2 Sample(Vector2 tangent)
    {
        //float angle = Mathf.Atan2(Mathf.Abs(tangent.y), Mathf.Abs(tangent.x));
        float angle = Mathf.Atan(tangent.y / tangent.x);
        //Debug.Log(angle);
        return angle >= 0 ? new Vector2(1, 1) : new Vector2(-1, 1);
    }

    public VertJacobian SpatialDerivative()
    {
        return new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };
    }
}

public enum PotentialType
{
    VectorField
}