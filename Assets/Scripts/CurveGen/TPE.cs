using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public class TPE
{
    static TPE instance;
    public static TPE GetInstance { get { if (instance == null) instance = new TPE(); return instance; } }

    float alpha { get { return EnergyCurve.alpha; } }
    float beta { get { return EnergyCurve.beta; } }

    /// <summary>
    /// Computes the energy gradients of all Vertices compared to all other vertices and adds them to the matrix for each Vertex
    /// </summary>
    public void FillGradientVectorDirect(EnergyCurve curve, Matrix<float> gradients)
    {
        int numVerts = curve.NumVerts();
        for (int i = 0; i < numVerts; i++)
        {
            for (int j = 0; j < numVerts; j++)
            {
                FillGradientSingle(curve, gradients, i, j);
            }
        }
    }

    public void FillGradientVectorBH(EnergyCurve curve, BVHNode2D treeRoot, Matrix<float> output)
    {
        // The single energy term (i, j) affects six vertices:
        // (i_prev, i, i_next, j_prev, j, j_next)
        // We can restructure the computation as follows:
        // for each single 1-ring (i, i_prev, i_next), accumulate
        // the contributinos from the gradients of both termns (i, j) and (j, i)
        int numVerts = curve.NumVerts();
        Matrix<float> partialOutput = output;

        for (int i = 0; i < numVerts; i++)
        {
            CurveVertex i_pt = curve.verts[i];
            treeRoot.AccumulateTPEGradient(partialOutput, i_pt, curve, alpha, beta);
        }

        output += partialOutput;
    }

    /// <summary>
    /// Computes the gradients of the energy from two vertices. Adds these two gradients to the vertices' gradients.
    /// </summary>
    void FillGradientSingle(EnergyCurve curve, Matrix<float> gradients, int i, int j)
    {
        if (i == j) return;
        CurveVertex i_pt = curve.verts[i];
        CurveVertex j_pt = curve.verts[j];

        // Add i and neighbors of i
        List<CurveVertex> i_pts = new List<CurveVertex>();
        i_pts.Add(i_pt);
        for (int e = 0; e < i_pt.NumEdges(); e++)
            i_pts.Add(i_pt.Edge(e).Opposite(i_pt));

        // Add j and neighbors of j
        List<CurveVertex> j_pts = new List<CurveVertex>();
        j_pts.Add(j_pt);
        for (int e = 0; e < j_pt.NumEdges(); e++)
            j_pts.Add(j_pt.Edge(e).Opposite(j_pt));

        // Differentiate wrt neighbors of i
        // Exp: Basically we calculate the gradient of the energy between i and j. We do this for i and his neighbors.
        // This helps us understand how the energy between i and j would change if we were to move i.
        foreach (CurveVertex i_n in i_pts)
            CurveGenUtils.AddToRow(gradients, i_n.GlobalIndex(), TpeGrad(i_pt, j_pt, i_n));
        // Differentiate wrt neighbors of j
        foreach (CurveVertex j_n in j_pts)
        {
            bool noOverlap = true;
            // Only compute this derivative if j_n is not already included in one of the previous pairs
            foreach (CurveVertex i_n in i_pts)
            {
                if (i_n.Equals(j_n))
                    noOverlap = false;
            }
            if (noOverlap)
                CurveGenUtils.AddToRow(gradients, j_n.GlobalIndex(), TpeGrad(i_pt, j_pt, j_n));
        }
    }

    /// <summary>
    /// Computes the gradient of the TPE-Kernel (K_f(x, y) dy dy) with respect to the position of the vertex "wrt"
    /// "wrt" just constrols whether the gradient is positive or negative.
    /// ToDo: Really understand what this method does in contrast to TpeGradKf()
    /// </summary>
    /// <param name="x"></param>
    /// <param name="y"></param>
    /// <param name="wrt"></param>
    /// <returns></returns>
    public Vector2 TpeGrad(CurveVertex x, CurveVertex y, CurveVertex wrt)
    {
        if (x.NumEdges() > 2) return Vector2.zero;

        // First, calculate the gradient of K_f(x, y) wrt "wrt"
        Vector2 grad_Kf = TpeGradKf(x, y, wrt);
        float Kf = TpeKf(x, y);

        float l_x = x.AvgLength();
        float l_y = y.AvgLength();

        // d/dy of area(x)
        Vector2 grad_lx = LengthWrtVert(x, wrt);
        // d/dy of area(y)
        Vector2 grad_ly = LengthWrtVert(y, wrt);
        // Evaluate the product rule for dx*dy
        Vector2 prod_rule = grad_lx * l_y + l_x * grad_ly;

        // Evaluate the product rule for k dx dy
        return grad_Kf * l_x * l_y + Kf * prod_rule;
    }

    public Vector2 TpeGrad(TangentMassPoint x, CurveVertex y, CurveVertex wrt)
    {
        // First, calculate the gradient of K_f(x, y)
        Vector2 grad_Kf = TpeGradKf(x, y, wrt);
        float Kf = TpeKfPts(x.point, y.Position(), x.tangent);

        float l_x = x.mass;
        float l_y = y.AvgLength();

        // Area gradient for x depends on whether the mass point is distant
        Vector2 grad_lx = Vector2.zero;
        if (x.GetTPEPointType() == TPEPointType.Point)
            grad_lx = LengthWrtVert(x.curvePt, wrt);
        else if (x.GetTPEPointType() == TPEPointType.Edge)
        {
            if (wrt.Equals(x.curvePt))
            {
                grad_lx = (x.curvePt.Position() - x.curvePt2.Position());
                grad_lx = grad_lx.normalized;
            }
            else if (wrt.Equals(x.curvePt2))
            {
                grad_lx = (x.curvePt2.Position() - x.curvePt.Position());
                grad_lx = grad_lx.normalized;
            }
        }
        // d/dy of area(y)
        Vector2 grad_ly = LengthWrtVert(y, wrt);
        // Evaluate the product rule for dy*dy
        Vector2 prod_rule = grad_lx * l_y + l_x * grad_ly;

        // Evalue the product rule for k dx dy
        return grad_Kf * l_x * l_y + Kf * prod_rule;
    }

    public Vector2 TpeGrad(CurveVertex x, TangentMassPoint y, CurveVertex wrt)
    {
        // First, calculate the gradient of K_f(x, y)
        Vector2 grad_Kf = TpeGradKf(x, y, wrt);
        float Kf = TpeKfPts(x.Position(), y.point, x.Tangent());

        float l_x = x.AvgLength();
        float l_y = y.mass;

        // Area gradient for x depends on whether the mass point is distant
        Vector2 grad_ly = Vector2.zero;
        if (y.GetTPEPointType() == TPEPointType.Point)
            grad_ly = LengthWrtVert(y.curvePt, wrt);
        else if (y.GetTPEPointType() == TPEPointType.Edge)
        {
            if (wrt.Equals(y.curvePt))
            {
                grad_ly = (y.curvePt.Position() - y.curvePt2.Position());
                grad_ly = grad_ly.normalized;
            }
            else if (wrt.Equals(y.curvePt2))
            {
                grad_ly = (y.curvePt2.Position() - y.curvePt.Position());
                grad_ly = grad_ly.normalized;
            }
        }
        // d/dy of area(y)
        Vector2 grad_lx = LengthWrtVert(x, wrt);
        // Evaluate the product rule for dy*dy
        Vector2 prod_rule = grad_lx * l_y + l_x * grad_ly;

        // Evalue the product rule for k dx dy
        return grad_Kf * l_x * l_y + Kf * prod_rule;
    }

    /// <summary>
    /// Tangent-Point-Energy Gradient Kernel-function
    /// </summary>
    Vector2 TpeGradKf(CurveVertex i, CurveVertex j, CurveVertex wrt)
    {
        if (i.Equals(j)) return Vector2.zero;

        // Get positions and displacement vectors
        Vector2 disp = i.Position() - j.Position();
        Vector2 T_i = i.Tangent();
        Vector2 unit_disp = disp.normalized; // Normalized displacement direction

        // Evaulate projection onto normal plane, v - <v,T> * T
        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float numer = Mathf.Pow(normal_proj.magnitude, alpha); // Numerator of energy is norm of projection ^ alpha
        float denom = Mathf.Pow(disp.magnitude, beta); // Denominator of energy is distance between points ^ beta
        // Our Kernel function (Eq. 3) would now be numer / denom. However, we want its derivative.

        // Derivative of numerator
        Vector2 deriv_numer = GradNormProjAlpha(i, j, wrt);

        // Derivative of denominator
        Vector2 deriv_denom = Vector2.zero;
        if (wrt.Equals(i))
            deriv_denom = beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;
        else if (wrt.Equals(j))
            deriv_denom = -beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;

        // Quotient rule for A / B
        Vector2 total = (deriv_numer * denom - deriv_denom * numer) / (denom * denom);
        return total;
    }

    Vector2 TpeGradKf(TangentMassPoint i, CurveVertex j, CurveVertex wrt)
    {
        Vector2 disp = i.point - j.Position();
        Vector2 T_i = i.tangent;
        Vector2 unit_disp = disp.normalized;

        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float numer = Mathf.Pow(normal_proj.magnitude, alpha);
        float denom = Mathf.Pow(disp.magnitude, beta);

        Vector2 deriv_numer = GradNormProjAlpha(i, j, wrt);

        Vector2 deriv_denom = Vector2.zero;
        if (i.GetTPEPointType() == TPEPointType.Point && wrt.Equals(i.curvePt))
            deriv_denom = beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;
        else if (i.GetTPEPointType() == TPEPointType.Edge && (wrt.Equals(i.curvePt) || wrt.Equals(i.curvePt2)))
            deriv_denom = (beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp) * 0.5f;
        if (wrt.Equals(j))
            deriv_denom += -beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;

        return (deriv_numer * denom - deriv_denom * numer) / (denom * denom);
    }

    Vector2 TpeGradKf(CurveVertex i, TangentMassPoint j, CurveVertex wrt)
    {
        Vector2 disp = i.Position() - j.point;
        Vector2 T_i = i.Tangent();
        Vector2 unit_disp = disp.normalized;

        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float numer = Mathf.Pow(normal_proj.magnitude, alpha);
        float denom = Mathf.Pow(disp.magnitude, beta);

        Vector2 deriv_numer = GradNormProjAlpha(i, j, wrt);

        Vector2 deriv_denom = Vector2.zero;
        if (wrt.Equals(i))
            deriv_denom = beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;
        if (j.GetTPEPointType() == TPEPointType.Point && wrt.Equals(j.curvePt))
            deriv_denom += -beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp;
        else if (j.GetTPEPointType() == TPEPointType.Edge && (wrt.Equals(j.curvePt) || wrt.Equals(j.curvePt2)))
            deriv_denom += (-beta * Mathf.Pow(disp.magnitude, beta - 1) * unit_disp) * 0.5f;

        return (deriv_numer * denom - deriv_denom * numer) / (denom * denom);
    }

    /// <summary>
    /// Kernel-function for the TPE, Equation (3). Returns the energy generated by i and j. 
    float TpeKf(CurveVertex i, CurveVertex j)
    {
        if (i.Equals(j) || i.NumEdges() > 2) return 0;

        Vector2 disp = i.Position() - j.Position();
        Vector2 t_i = i.Tangent();

        Vector2 normal_proj = disp - Vector2.Dot(disp, t_i) * t_i;
        float numer = Mathf.Pow(normal_proj.magnitude, alpha);
        float denom = Mathf.Pow(disp.magnitude, beta);
        return numer / denom;
    }

    private float TpeKfPts(Vector2 p_x, Vector2 p_y, Vector2 tangent_x)
    {
        Vector2 disp = p_x - p_y;
        Vector2 n_proj = disp - Vector2.Dot(disp, tangent_x) * tangent_x;
        float numer = Mathf.Pow(n_proj.magnitude, alpha);
        float denom = Mathf.Pow(disp.magnitude, beta);
        return numer / denom;
    }

    public float TpePair(CurveVertex i, CurveVertex j)
    {
        float kfxy = TpeKf(i, j);
        float l_x = i.AvgLength();
        float l_y = i.AvgLength();
        return l_x * l_y * kfxy;
    }

    public float TpePairPts(Vector2 p_x, Vector2 p_y, Vector2 tangent_x, float l_x, float l_y)
    {
        float kfxy = TpeKfPts(p_x, p_y, tangent_x);
        return l_x * l_y * kfxy;
    }

    public float TpeTotal(EnergyCurve curve)
    {
        int nVerts = curve.NumVerts();
        float sumEnergy = 0;

        for (int i = 0; i < nVerts; i++)
        {
            for (int j = 0; j < nVerts; j++)
            {
                CurveVertex pt_i = curve.verts[i];
                CurveVertex pt_j = curve.verts[j];
                sumEnergy += TpePair(pt_i, pt_j);
            }
        }

        return sumEnergy;
    }

    /// <summary>
    /// Returns the "length" of one Vert with respect to another vert "wrt".
    /// If lengthVert != wrt, this is just (lengthVert.pos - wrt.pos).normalized * 0.5f
    /// Otherwise it's all the outgoing directions * 0.5f
    /// </summary>
    /// <param name="lengthVert"></param>
    /// <param name="wrt"></param>
    /// <returns></returns>
    Vector2 LengthWrtVert(CurveVertex lengthVert, CurveVertex wrt)
    {
        // If differentiating wrt self, need to consider both sides
        if (lengthVert == wrt)
        {
            Vector2 sumDirections = Vector2.zero;
            Vector2 center = lengthVert.Position();
            for (int e = 0; e < lengthVert.NumEdges(); e++)
            {
                Vector2 other = lengthVert.Edge(e).Opposite(lengthVert).Position();
                Vector2 outward = other - center;
                outward = outward.normalized;
                sumDirections += outward;
            }
            // Gradient is half the sum of unit vectors along outgoing edges.
            return -0.5f * sumDirections;
        }
        // Otherwise, only consider the one edge between other and vert
        else
        {
            for (int e = 0; e < lengthVert.NumEdges(); e++)
            {
                CurveVertex other = lengthVert.Edge(e).Opposite(lengthVert);
                if (other == wrt)
                {
                    Vector2 outward = other.Position() - lengthVert.Position();
                    outward = outward.normalized;
                    return 0.5f * outward;
                }
            }
            return Vector2.zero;
        }
    }

    /// <summary>
    /// Basically differentiates the numerater of the Kernel-Function (Eq. 3)
    /// </summary>
    /// <param name="i"></param>
    /// <param name="j"></param>
    /// <param name="wrt"></param>
    /// <returns></returns>
    Vector2 GradNormProjAlpha(CurveVertex i, CurveVertex j, CurveVertex wrt)
    {
        Vector2 disp = i.Position() - j.Position();
        Vector2 T_i = i.Tangent();
        // Projection
        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float proj_len = normal_proj.magnitude; // same as T_i x disp = CurveGenUtils.Cross(T_i, disp)

        // If the displacement is usually exactly perpendicular to the tangent, then the contribution is exactly 0;
        if (proj_len < 1e-10) return Vector2.zero;

        // Derivative of |f(x) - ...|^alpha = alpha * |f(x) - ...|^(alpha-1)
        float alpha_deriv = alpha * Mathf.Pow(proj_len, alpha - 1);
        // Normalized vector of projection onto normal vector
        Vector2 proj_normalized = normal_proj / proj_len;

        VertJacobian deriv_disp = new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };
        if (wrt.Equals(i))
        {
            deriv_disp.directional_x = new Vector2(1, 0);
            deriv_disp.directional_y = new Vector2(0, 1);
        }
        else if (wrt.Equals(j))
        {
            deriv_disp.directional_x = new Vector2(-1, 0);
            deriv_disp.directional_y = new Vector2(0, -1);
        }

        // Derivative of <f(x) - f(y), T> * T
        VertJacobian deriv_T_inner = GradTangentProj(i, j, wrt); // deriv(<d,T>*T)
        VertJacobian deriv_N_Proj = deriv_disp - deriv_T_inner; // deriv(disp) - deriv(<d,T>*T)

        //Debug.Log(proj_len);
        //float cross = CurveGenUtils.Cross(T_i, disp);
        //Debug.Log(cross);

        Vector2 total = alpha_deriv * deriv_N_Proj.LeftMultiply(proj_normalized);
        return total;
    }

    Vector2 GradNormProjAlpha(TangentMassPoint i, CurveVertex j, CurveVertex wrt)
    {
        Vector2 disp = i.point - j.Position();
        Vector2 T_i = i.tangent;
        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float proj_len = normal_proj.magnitude;

        if (proj_len < 1e-10) return Vector2.zero;

        float alpha_deriv = alpha * Mathf.Pow(proj_len, alpha - 1);
        Vector2 proj_normalized = normal_proj / proj_len;

        VertJacobian deriv_disp = new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };
        if (i.GetTPEPointType() == TPEPointType.Point && wrt.Equals(i))
        {
            deriv_disp.directional_x = new Vector2(1, 0);
            deriv_disp.directional_y = new Vector2(0, 1);
        }
        if (i.GetTPEPointType() == TPEPointType.Edge && (wrt.Equals(i.curvePt) || wrt.Equals(i.curvePt2)))
        {
            deriv_disp.directional_x = new Vector2(0.5f, 0);
            deriv_disp.directional_y = new Vector2(0, 0.5f);
        }
        else if (wrt.Equals(j))
        {
            deriv_disp.directional_x = new Vector2(-1, 0);
            deriv_disp.directional_y = new Vector2(0, -1);
        }

        VertJacobian deriv_T_inner = GradTangentProj(i, j, wrt);
        VertJacobian deriv_N_Proj = deriv_disp - deriv_T_inner;

        return alpha_deriv * deriv_N_Proj.LeftMultiply(proj_normalized);
    }

    Vector2 GradNormProjAlpha(CurveVertex i, TangentMassPoint j, CurveVertex wrt)
    {
        Vector2 disp = i.Position() - j.point;
        Vector2 T_i = i.Tangent();
        Vector2 normal_proj = disp - Vector2.Dot(disp, T_i) * T_i;
        float proj_len = normal_proj.magnitude;

        if (proj_len < 1e-10) return Vector2.zero;

        float alpha_deriv = alpha * Mathf.Pow(proj_len, alpha - 1);
        Vector2 proj_normalized = normal_proj / proj_len;

        VertJacobian deriv_disp = new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };
        if (wrt.Equals(i))
        {
            deriv_disp.directional_x = new Vector2(1, 0);
            deriv_disp.directional_y = new Vector2(0, 1);
        }
        else if (j.GetTPEPointType() == TPEPointType.Point && wrt.Equals(j))
        {
            deriv_disp.directional_x = new Vector2(-1, 0);
            deriv_disp.directional_y = new Vector2(0, -1);
        }
        else if (j.GetTPEPointType() == TPEPointType.Edge && (wrt.Equals(j.curvePt) || wrt.Equals(j.curvePt2)))
        {
            deriv_disp.directional_x = new Vector2(-0.5f, 0);
            deriv_disp.directional_y = new Vector2(0, -0.5f);
        }

        VertJacobian deriv_T_inner = GradTangentProj(i, j, wrt);
        VertJacobian deriv_N_Proj = deriv_disp - deriv_T_inner;

        return alpha_deriv * deriv_N_Proj.LeftMultiply(proj_normalized);
    }

    /// <summary>
    /// // ddx(<d,T>*T)
    /// </summary>
    VertJacobian GradTangentProj(CurveVertex i, CurveVertex j, CurveVertex wrt)
    {
        // Differentiate the inner product
        Vector2 disp = i.Position() - j.Position();
        Vector2 T_i = i.Tangent();
        float disp_dot_T = Vector2.Dot(disp, T_i);

        Vector2 inner_deriv_A_B = Vector2.zero;                     // ddx(d) * T
        if (wrt.Equals(i))
            inner_deriv_A_B = T_i;
        else if (wrt.Equals(j))
            inner_deriv_A_B = -T_i;
        VertJacobian deriv_T = VertexTangentWrtVert(i, wrt);        // ddx(T)
        Vector2 inner_A_deriv_B = deriv_T.LeftMultiply(disp);       // d * ddx(T)
        Vector2 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;    // ddx(<d,T>) = ddx(d) * T + d * ddx(T)

        // Now use product rule for <f(x) - f(y), T> * T
        VertJacobian deriv_A_B = VertJacobian.OuterProductToJacobian(T_i, deriv_inner); // ddx(<d,T>) * T
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;                                  // <d,T> * ddx(T)

        return deriv_A_B + A_deriv_B;   // deriv(<d,T>) * T + <d,T> * deriv(T)
    }

    VertJacobian GradTangentProj(TangentMassPoint i, CurveVertex j, CurveVertex wrt)
    {
        Vector2 disp = i.point - j.Position();
        Vector2 T_i = i.tangent;
        float disp_dot_T = Vector2.Dot(disp, T_i);

        Vector2 inner_deriv_A_B = Vector2.zero;
        if (wrt.Equals(j))
            inner_deriv_A_B = -T_i;
        VertJacobian deriv_T = VertexTangentWrtVert(i.curvePt, wrt);
        Vector2 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
        Vector2 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

        VertJacobian deriv_A_B = VertJacobian.OuterProductToJacobian(T_i, deriv_inner);
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;

        return deriv_A_B + A_deriv_B;
    }

    VertJacobian GradTangentProj(CurveVertex i, TangentMassPoint j, CurveVertex wrt)
    {
        Vector2 disp = i.Position() - j.point;
        Vector2 T_i = i.Tangent();
        float disp_dot_T = Vector2.Dot(disp, T_i);

        Vector2 inner_deriv_A_B = Vector2.zero;
        if (wrt.Equals(i))
            inner_deriv_A_B = T_i;
        VertJacobian deriv_T = VertexTangentWrtVert(i, wrt);
        Vector2 inner_A_deriv_B = deriv_T.LeftMultiply(disp);
        Vector2 deriv_inner = inner_deriv_A_B + inner_A_deriv_B;

        VertJacobian deriv_A_B = VertJacobian.OuterProductToJacobian(T_i, deriv_inner);
        VertJacobian A_deriv_B = disp_dot_T * deriv_T;

        return deriv_A_B + A_deriv_B;
    }

    /// <summary>
    /// ddx(T) where T = (T1 + T2) / |T1 + T2|. Derivative of Vertex-Tangent with respect to Vert
    /// </summary>
    VertJacobian VertexTangentWrtVert(CurveVertex tangentVert, CurveVertex wrtVert)
    {
        if (tangentVert == null || wrtVert == null)
            return new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };

        if (tangentVert.NumEdges() != 2)
        {
            if (tangentVert.NumEdges() == 1)
            {
                CurveEdge edge = tangentVert.Edge(0);
                Vector2 tangent = edge.Tangent();
                // Derivative of T
                return EdgeTangentWrtVert(edge, wrtVert);
            }
            else
            {
                return new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };
            }
        }

        CurveEdge prevEdge = tangentVert.Edge(0);
        CurveEdge nextEdge = tangentVert.Edge(1);

        Vector2 prevTangent = prevEdge.Tangent();
        Vector2 nextTangent = nextEdge.Tangent();

        Vector2 sumTangents = prevTangent + nextTangent;    //  T1+T2
        float normSum = sumTangents.magnitude;              // |T1+T2|

        VertJacobian derivSumTs = EdgeTangentWrtVert(prevEdge, wrtVert) + EdgeTangentWrtVert(nextEdge, wrtVert);    // ddx(T1+T2) = ddx(T1)+ddx(T2)

        Vector2 derivNorm = derivSumTs.LeftMultiply(tangentVert.Tangent());         // ddx(|T1+T2|) = ddx(T1+T2) * T
        VertJacobian deriv_A_B = derivSumTs * normSum;                                          // ddx(T1+T2) * |T1+T2|
        VertJacobian A_deriv_B = VertJacobian.OuterProductToJacobian(sumTangents, derivNorm);   // (T1+T2) * ddx(|T1+T2|)

        return (deriv_A_B - A_deriv_B) / (normSum * normSum); // Quotient rule on (T1+T2) / |T1+T2|
    }

    /// <summary>
    /// ddx(T). Derviative of Edge-Tangent with respect to Vert
    /// </summary>
    public VertJacobian EdgeTangentWrtVert(CurveEdge edge, CurveVertex wrtVert)
    {
        // Get positions
        CurveVertex prevVert = edge.GetPrevVertex();
        CurveVertex nextVert = edge.GetNextVertex();
        Vector2 v_h = prevVert.Position();
        Vector2 v_i = nextVert.Position();

        if (!wrtVert.Equals(prevVert) && !wrtVert.Equals(nextVert)) // Logic doesn't work if wrt is not part of the edge
            return new VertJacobian { directional_x = Vector2.zero, directional_y = Vector2.zero };

        Vector2 v_tangent = v_i - v_h; // not normalized
        float v_norm = v_tangent.magnitude;
        Vector2 v_normalized = v_tangent / v_norm; // ddx(T)

        VertJacobian I = new VertJacobian { directional_x = new Vector2(1, 0), directional_y = new Vector2(0, 1) };

        VertJacobian deriv_A_B = I * v_norm; // Scaling the Jacobian        
        VertJacobian A_deriv_B = VertJacobian.OuterProductToJacobian(v_tangent, v_normalized);
        VertJacobian deriv = (deriv_A_B - A_deriv_B) / (v_norm * v_norm); // Quotient rule

        // If we're differentiating the tail vertex, the derivative is negative
        if (wrtVert == prevVert) return -1 * deriv;
        else return deriv;
    }

    public float TPE_Kf_pts_sym(Vector2 p_x, Vector2 p_y, Vector2 tangent_x, Vector2 tangent_y,
        float alpha, float beta) // alpha and beta may have other values than usual
    {
        Vector2 disp = p_x - p_y;
        Vector2 n_proj_x = disp - Vector2.Dot(disp, tangent_x) * tangent_x;
        Vector2 n_proj_y = disp - Vector2.Dot(disp, tangent_y) * tangent_y;

        float numer_x = Mathf.Pow(n_proj_x.magnitude, alpha);
        float numer_y = Mathf.Pow(n_proj_y.magnitude, alpha);
        float denom = Mathf.Pow(disp.magnitude, beta);

        return 0.5f * (numer_x + numer_y) / denom;
    }

    internal Vector2 EdgeLengthWrtVert(CurveEdge edge, CurveVertex wrt)
    {
        if (!edge.GetPrevVertex().Equals(wrt) && !edge.GetNextVertex().Equals(wrt))
            return Vector2.zero;

        Vector2 wrt_pos = wrt.Position();
        Vector2 opp_pos = edge.Opposite(wrt).Position();

        // To increase the edge length, we want to move away from the opposite vertex
        Vector2 outward = wrt_pos - opp_pos;
        return outward.normalized;
    }
}
