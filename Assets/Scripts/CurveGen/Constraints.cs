using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class VariableConstraintSet
{
    readonly EnergyCurve curve;

    public List<Constraint> constraints;
    public List<ConstraintType> ConstraintTypes { get { return constraints.Select(c => c.type).ToList(); } }

    #region Init

    public VariableConstraintSet(EnergyCurve curve, List<ConstraintType> types)
    {
        this.curve = curve;

        constraints = new List<Constraint>();
        foreach (ConstraintType type in types)
        {
            switch (type)
            {
                case ConstraintType.Length: constraints.Add(new LengthConstraint(curve)); break;
                case ConstraintType.Barycenter: constraints.Add(new BarycenterConstraint(curve)); break;
                default: throw new ArgumentException("Supplied type is not implemented yet");
            }
        }
    }

    #endregion

    public void AddTriplets(List<Triplet> triplets)
    {
        int startIndex = 0;
        foreach (Constraint constraint in constraints)
        {
            constraint.AddTriplets(triplets, startIndex);
            startIndex += NumRowsForConstraint(constraint.type);
        }
    }

    public void SetTargetValues(Vector<float> targets)
    {
        int startIndex = 0;
        foreach (Constraint constraint in constraints)
        {
            constraint.SetTargets(targets, startIndex);
            startIndex += NumRowsForConstraint(constraint.type);
        }
    }

    /// <summary>
    /// Given a vector b with NumConstraintRows() entries, and given the target
    /// values of the constraint function, fills b with the corresponding negative
    /// values of the function.
    /// </summary>
    public void NegativeConstraintValues(Vector<float> b, Vector<float> targets)
    {
        int startIndex = 0;
        foreach (Constraint constraint in constraints)
        {
            constraint.NegativeViolation(b, targets, startIndex);
            startIndex += NumRowsForConstraint(constraint.type);
        }
    }

    #region Helper

    public int NumConstraintRows()
    {
        int rows = 0;
        foreach (Constraint constraint in constraints)
        {
            rows += NumRowsForConstraint(constraint.type);
        }
        return rows;
    }

    public int NumExpectedCols()
    {
        return curve.NumVerts() * 2;
    }


    public int StartIndexOfConstraint(ConstraintType type)
    {
        int startIndex = 0;
        foreach (Constraint constraint in constraints)
        {
            if (type == constraint.type) return startIndex;
            else startIndex += NumRowsForConstraint(constraint.type);
        }
        return -1;
    }

    public int NumRowsForConstraint(ConstraintType type)
    {
        return type switch
        {
            ConstraintType.Barycenter => 2,
            ConstraintType.Length => curve.NumEdges(),
            _ => throw new ArgumentException("Called NumRowsForConstraint on an unimplemented constraint type."),
        };
    }

    /// <summary>
    /// Given a dense matrix sps1 (first matrix of saddle point system) of size SaddleNumRows() x SaddleNumRows()
    /// (which is the size of a full saddle matrix), fills the blocks
    /// for C and C^T at the bottom and right edges of the matrix.
    /// The input matrix must already have this size.
    /// </summary>
    public void FillDenseBlock(Matrix<float> sps1)
    {
        List<Triplet> triplets = new();
        AddTriplets(triplets);
        int offset = NumExpectedCols();
        foreach (Triplet t in triplets)
        {
            // Copy into lower-left block
            //Debug.Log("Triplet: (" + t.row + ", " + t.column + ") -> " + t.value);
            sps1[offset + t.row, t.column] = t.value;
            sps1[t.column, offset + t.row] = t.value;
        }
    }

    /// <summary>
    /// Fills the given vector b with the negated current values of this constraint function.
    /// Entries are written starting from the index of the argument "offfset".
    /// </summary>
    internal float FillConstraintValues(Vector<float> b, Vector<float> targets, int offset)
    {
        // Fill values of constraint function in a generic way
        int numConstraintRows = NumConstraintRows();
        Vector<float> b_constraints = Vector<float>.Build.Dense(numConstraintRows);
        NegativeConstraintValues(b_constraints, targets);
        b.SetSubVector(offset, numConstraintRows, b_constraints);

        return (float)b_constraints.InfinityNorm(); // Summe der absoluten Werte
    }

    /// <summary>
    /// Given a vector containing the target values for the constraint function,
    /// this function updates each target value with the current corresponding target value
    /// for the constraint.
    /// </summary>
    internal void UpdateTargetValues(Vector<float> targets)
    {
        int numConstraints = NumConstraintRows();
        if (targets.Count != numConstraints)
            targets = Vector<float>.Build.Dense(numConstraints);
        SetTargetValues(targets);
    }

    #endregion
}

[Serializable]
public enum ConstraintType
{
    Barycenter, Length
}

public struct Triplet
{
    public int row, column;
    public float value;
}


public abstract class Constraint
{
    public ConstraintType type;
    public EnergyCurve curve;

    public Constraint(EnergyCurve curve, ConstraintType type)
    {
        this.curve = curve;
        this.type = type;
    }

    public abstract void SetTargets(Vector<float> targets, int rowStart);

    public abstract void AddTriplets(List<Triplet> triplets, int rowStart);

    public abstract void NegativeViolation(Vector<float> b, Vector<float> targets, int rowStart);
}


public class LengthConstraint : Constraint
{
    public LengthConstraint(EnergyCurve curve) : base(curve, ConstraintType.Length)
    {
    }

    public override void AddTriplets(List<Triplet> triplets, int rowStart)
    {
        int numEdges = curve.NumEdges();

        // Add the edge length rows
        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge edge = curve.edges[i];
            CurveVertex pt1 = edge.GetPrevVertex();
            CurveVertex pt2 = edge.GetNextVertex();

            // This is the gradient of edge length wrt pt1;
            // the gradient wrt pt2 is just negative of this.
            Vector2 grad1 = pt1.Position() - pt2.Position();
            grad1 = grad1.normalized;

            int j1 = pt1.GlobalIndex();
            int j2 = pt2.GlobalIndex();
            int start = i + rowStart;

            // Write the two gradient entries for pt1 into the row
            triplets.Add(new Triplet { row = start, column = 2 * j1, value = grad1.x });
            triplets.Add(new Triplet { row = start, column = 2 * j1 + 1, value = grad1.y });

            // Similarly, write the two gradient entries for pt2 into the same row
            triplets.Add(new Triplet { row = start, column = 2 * j2, value = -grad1.x });
            triplets.Add(new Triplet { row = start, column = 2 * j2 + 1, value = -grad1.y });
        }
    }

    public override void NegativeViolation(Vector<float> b, Vector<float> targets, int rowStart)
    {
        int numEdges = curve.NumEdges();

        // For each edge, fill in its derivatives from target length
        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge e_i = curve.edges[i];
            float currLen = e_i.Length();
            int id = e_i.GetIndex();
            float negErr = targets[rowStart + id] - currLen;
            b[rowStart + id] = negErr;
        }
    }

    public override void SetTargets(Vector<float> targets, int rowStart)
    {
        int numEdges = curve.NumEdges();
        for (int i = 0; i < numEdges; i++)
            targets[rowStart + i] = curve.edges[i].Length();
    }
}


public class BarycenterConstraint : Constraint
{
    public BarycenterConstraint(EnergyCurve curve) : base(curve, ConstraintType.Barycenter)
    {
    }

    public override void AddTriplets(List<Triplet> triplets, int rowStart)
    {
        int numVerts = curve.NumVerts();
        float totalLength = curve.TotalLength();

        // Fill a single row with normalized vertex weights
        for (int i = 0; i < numVerts; i++)
        {
            float wt = curve.verts[i].AvgLength() / totalLength;
            triplets.Add(new Triplet { row = rowStart + 0, column = 2 * i, value = wt });
            triplets.Add(new Triplet { row = rowStart + 1, column = 2 * i + 1, value = wt });
        }
    }

    public override void NegativeViolation(Vector<float> b, Vector<float> targets, int rowStart)
    {
        Vector2 barycenter = curve.Barycenter();

        // Fill in difference from barycenter to target points
        b[rowStart + 0] = targets[rowStart + 0] - barycenter.x;
        b[rowStart + 1] = targets[rowStart + 1] - barycenter.y;
    }

    public override void SetTargets(Vector<float> targets, int rowStart)
    {
        Vector2 barycenter = curve.Barycenter();

        // Set the two barycenter target entries to 0
        targets[rowStart + 0] = barycenter.x;
        targets[rowStart + 1] = barycenter.y;
    }
}