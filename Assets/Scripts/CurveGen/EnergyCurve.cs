using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using Random = UnityEngine.Random;

public class EnergyCurve : Curve
{
    #region Variables

    public Spline spline = null;
    Matrix<float> originalPositions;

    public List<Obstacle> obstacles; // ToDo: Add obstacles as type of potential
    public List<Potential> potentials;

    public VariableConstraintSet constraintSet;
    Vector<float> constraintTargets;

    public const float alpha = 3f, beta = 6f;

    float lastStepSize;
    /// <summary>
    /// Smallest possible step size
    /// </summary>
    const float lsStepThreshold = 1e-15f;
    /// <summary>
    /// If the "length" (Frobeniusnorm) of the gradient-matrix is less than a certain threshold (we're near a minimum), this parameter will be true
    /// </summary>
    public bool normZero;
    const float backproj_threshold = 1e-4f;
    float startTime;

    // Lengths
    float lengthScaleStep;
    public float initialAvgLength;
    float targetLength;
    bool scalingLength;
    bool deacObsAfterScaling;
    bool rotateAfterScaling;
    bool noRepulsionAfterScaling;

    // For Reproducability (Reset)
    List<Vector2> origPosList;
    readonly float innerObstacleRadius = 1;

    #region Serialization
    List<S_ObstacleConfig> s_obstacles;
    #endregion

    bool Repulsion { get { return scalingLength || !noRepulsionAfterScaling; } }

    #endregion

    #region Inits

    /// <summary>
    /// Generate the curve. Called by Editor.
    /// If genMode == Full, then use construct the curve as a whole and jus run LS.
    /// If genMode == Grow, construct a very simple curve that's going to be grown afterwards.
    /// </summary>
    public EnergyCurve(EnergyCurve_EditorConfig config)
    {
        curveClosed = config.curveClosed;

        #region Curve

        List<Vector2> posList = new();

        if (config.genModeConfig is GenModeConfig_Bezier bezierConfig)
        {
            List<Vector2> splinePoints = new();
            for (int i = 0; i < bezierConfig.numPoints; i++)
            {
                float alpha = Random.Range(-Mathf.PI, Mathf.PI);
                splinePoints.Add(bezierConfig.radius * Random.value * new Vector2(Mathf.Cos(alpha), Mathf.Sin(alpha)));
            }

            //points = FPDS.Sampling(Vector2.one * -radius, Vector2.one * radius, radius / 10, numPoints);
            spline = new(splinePoints)
            {
                AutoSetControlPoints = true
            };

            if (curveClosed) spline.ToggleClosed();

            posList = spline.CalculateEvenlySpacedPoints(bezierConfig.spacing, 10f);
        }
        else if (config.genModeConfig is GenModeConfig_Circular circularConfig)
        {
            for (int i = 0; i < circularConfig.numPoints; i++)
            {
                float alpha = (i / (float)circularConfig.numPoints) * Mathf.PI * 2f;
                posList.Add(circularConfig.radius * new Vector2(Mathf.Cos(alpha), Mathf.Sin(alpha)));
            }
            //points.ForEach(x => Debug.Log(x));
        }

        InitVertsEdgesFromPositions(posList);
        origPosList = new(posList);

        #endregion

        lastStepSize = 0;
        normZero = false;

        initialAvgLength = TotalLength() / NumEdges();
        targetLength = config.targetLengthScale * TotalLength();
        InitConstraints(config.constraints);

        #region Obstacles

        obstacles = new();
        if (config.obstacles != null)
        {
            foreach (ObstacleConfig obsConf in config.obstacles)
            {
                obstacles.Add(new(p_exp: beta - alpha, weight: obsConf.weight, numPoints: obsConf.numPoints, radius: obsConf.radius, center: obsConf.center));
            }
        }
        else
        {
            if (config.genModeConfig is GenModeConfig_Circular circularConfig)
            {
                s_obstacles = new();

                float innerRadius = circularConfig.radius;
                float outerRadius = innerRadius * 4;

                s_obstacles.Add(new S_ObstacleConfig { p_exp = beta - alpha, weight = 1, numPoints = (int)outerRadius * 3, radius = outerRadius, center = Vector2.zero });

                int numInnerObstacles = config.numObstacles; // - 1;

                for (int i = 0; i < numInnerObstacles; i++)
                {
                    float angle = i * 2 * Mathf.PI / numInnerObstacles;
                    float radius = Random.Range(innerRadius + innerObstacleRadius + 1, outerRadius);
                    Vector2 pos = new(Mathf.Cos(angle) * radius, Mathf.Sin(angle) * radius);

                    s_obstacles.Add(new S_ObstacleConfig { p_exp = beta - alpha, weight = 0.25f, numPoints = 10, radius = innerObstacleRadius, center = pos });
                }

                S_TryInitObstacles();
            }
            else throw new Exception("Algorithmically generating obstacles is only supported with Circular GenMode!");
        }

        this.deacObsAfterScaling = config.deacObsAfterScaling;
        this.rotateAfterScaling = config.rotateAfterScaling;
        this.noRepulsionAfterScaling = config.noRepulsionAfterScaling;

        #endregion

        InitPotentials(config.potentials);
    }

    private void InitPotentials(List<PotentialConfig> potentialConfigs)
    {
        potentials = new();
        foreach (PotentialConfig pConf in potentialConfigs)
        {
            switch (pConf.type)
            {
                case PotentialType.VectorField:
                    potentials.Add(new VectorFieldPotential(this, pConf.weight, pConf.delayed)); break;
                default: throw new ArgumentException("Passed PotentialType '" + pConf.type + "' is not implemented yet.");
            }
        }
    }

    public void InitConstraints(List<ConstraintType> constraintTypes)
    {
        constraintSet = new(this, constraintTypes);


        if (constraintTypes.Contains(ConstraintType.Length))
        {
            float avgLength = TotalLength() / NumEdges();
            lengthScaleStep = avgLength * 0.01f;
            scalingLength = true;
        }

        if (constraintTypes.Contains(ConstraintType.Barycenter))
        {
        }

        constraintTargets = Vector<float>.Build.Dense(constraintSet.NumConstraintRows());

        constraintSet.UpdateTargetValues(constraintTargets);
    }

    public void Reset(EnergyCurve_EditorConfig config)
    {
        InitVertsEdgesFromPositions(origPosList);

        lastStepSize = 0;
        normZero = false;
        targetLength = config.targetLengthScale * TotalLength();

        InitConstraints(config.constraints);
        InitPotentials(config.potentials);

        obstacles.ForEach(x => x.Enable());

        this.deacObsAfterScaling = config.deacObsAfterScaling;
        this.rotateAfterScaling = config.rotateAfterScaling;
        this.noRepulsionAfterScaling = config.noRepulsionAfterScaling;
    }

    #endregion

    #region Gradient Flow

    #region Measuring Time

    float GetTime() => Time.realtimeSinceStartup;
    float SetStartTime() => startTime = GetTime();
    string DurationString() => "(time: " + (GetTime() - startTime) + "s)";

    #endregion

    public Tuple<bool, bool> ComputeLineSearchStep(RepulsionType repulsionType, bool useBarnesHut, bool useBackproj)
    {
        return repulsionType switch
        {
            RepulsionType.Normal => ComputeLineSearchStep_Normal(useBarnesHut, useBackproj),
            RepulsionType.Sobolev => ComputeLineSearchStep_Sobolev(useBarnesHut, useBackproj),
            _ => throw new ArgumentException("Called ComputeLineSearchStep with wrong RepulsionType"),
        };
    }

    /// <summary>
    /// Compute one Iteration. Main Steps:
    /// (i) Assemble the gradient, (ii) apply constraints, (iii) do a line search step
    /// </summary>
    public Tuple<bool, bool> ComputeLineSearchStep_Normal(bool useBarnesHut, bool useBackproj)
    {
        int numVerts = NumVerts();
        Matrix<float> gradients = Matrix<float>.Build.Dense(numVerts, 2);

        // If applicable (if genMode == Grow), move constraint targets
        //SetStartTime();
        if (scalingLength) MoveLengthTowardsTarget();
        //Debug.Log("Moved Length towards target " + DurationString());

        #region Gradient Assembling
        // Assemble gradients, either exaclty or with Barnes-Hut.
        BVHNode2D treeRoot = null;
        if (useBarnesHut)
        {
            //SetStartTime();
            treeRoot = TPEBVH.GetInstance.CreateBVHFromCurve(this);
            //Debug.Log("Computed BH-Tree " + DurationString());
        }
        //SetStartTime();
        AddAllGradients(treeRoot, gradients);
        //Debug.Log("Assembled Gradients " + DurationString());
        //Debug.Log("GRADIENTS: " + gradients);
        #endregion

        #region Constraints
        if (constraintSet.constraints.Count > 0)
        {
            // Set up saddle matrix
            int numRows = constraintSet.NumConstraintRows() + constraintSet.NumExpectedCols();
            Matrix<float> A = Matrix<float>.Build.Dense(numRows, numRows);

            // Assemble constraint saddle matrix with identity in the upper-left corner
            Matrix<float> mass = Matrix<float>.Build.DenseIdentity(numVerts * 2, numVerts * 2);
            for (int i = 0; i < numVerts; i++)
            {
                float m = 1.0f / verts[i].AvgLength();
                mass[2 * i, 2 * i] = m;
                mass[2 * i + 1, 2 * i + 1] = m;
            }

            A.SetSubMatrix(0, 0, mass);

            constraintSet.FillDenseBlock(A);
            // Factorize it
            LU<float> lu = A.LU();

            //Debug.Log("Gradients before projection: " + gradients);

            // Project gradient onto constraint differential
            ProjectGradient(lu, gradients);

            //Debug.Log("Gradients after projection: " + gradients);

            //Debug.Log("NEW GRADIENTS: " + gradients);
        }
        #endregion

        #region Line Search
        //SetStartTime();
        float startLength = TotalLength();
        float gradNorm = (float)gradients.FrobeniusNorm();
        float stepSize = LineSearchStep(gradients, 1, treeRoot);
        //Debug.Log("Executed Line Search Step " + DurationString());
        //Debug.Log("New Length: " + startLength + " -> " + TotalLength());
        #endregion

        normZero = gradNorm < 1e-4f;
        lastStepSize = stepSize;
        return Tuple.Create((stepSize > lsStepThreshold), true);
    }

    /// <summary>
    /// (StepSobolevLS) Compute one Iteration using Sobolev. Main Steps:
    /// (i) Assemble Sobolev gradient and compute it, (ii) do one line search step with this gradient and
    /// (iii) correct for drift using backprojection
    /// </summary>
    public Tuple<bool, bool> ComputeLineSearchStep_Sobolev(bool useBarnesHut, bool useBackproj)
    {
        SetStartTime();

        int numVerts = NumVerts();
        Matrix<float> vertGradients = Matrix<float>.Build.Dense(numVerts, 2);

        //SetStartTime();

        // If applicable, move constraint targets
        if (scalingLength) MoveLengthTowardsTarget();
        if (scalingLength && TargetLengthReached())
        {
            Debug.Log("Target length reached. Not running length scaling");
            scalingLength = false;
            if (deacObsAfterScaling)
            {
                Debug.Log("Removing all obstacles");
                obstacles.ForEach(x => x.Disable());
            }

            bool shouldEnd = AfterScaling();
            //if (shouldEnd)
            //    return Tuple.Create(false, false);

        }

        Debug.Log("Moved Length towards target " + DurationString());

        #region Gradient Assembling
        // Assemble gradients, either exaclty or with Barnes-Hut.
        SetStartTime();
        BVHNode2D treeRoot = null;
        if (useBarnesHut)
        {
            treeRoot = TPEBVH.GetInstance.CreateBVHFromCurve(this);
            //Debug.Log("Computed BH-Tree " + DurationString());
        }
        AddAllGradients(treeRoot, vertGradients);
        Debug.Log("Assembled Gradients " + DurationString());
        Matrix<float> l2Gradients = vertGradients;
        #endregion

        #region Projection
        //SetStartTime();
        LU<float> lu = null;
        float soboDot = ProjectSoboSloboGradient(vertGradients, ref lu);
        //Debug.Log("Projected Gradient " + DurationString());
        if (soboDot == float.NaN)
        {
            Debug.Log("Sobolev projection produced NaN; aborting.");
            return Tuple.Create(false, true);
        }
        //else Debug.Log("Sobolev gradient norm: " + soboDot);
        #endregion

        float length1 = TotalLength();
        float energy1 = CurrentEnergy(treeRoot);

        float dot_acc = (float)(soboDot / (l2Gradients.FrobeniusNorm() * vertGradients.FrobeniusNorm()));
        //Debug.Log("dot_acc: " + dot_acc);

        #region Line Search
        // Take a line search step using this gradient
        //SetStartTime();

        bool resetStepSize = false;
        float stepSize = LineSearchStep(vertGradients, dot_acc, treeRoot, resetStepSize);
        //Debug.Log("Computed Line Search " + DurationString());
        #endregion

        if (UsingConstraint(ConstraintType.Length) && stepSize < lsStepThreshold)
            vertGradients = Matrix<float>.Build.Dense(numVerts, 2);

        #region Drift Correction with Backprojection
        if (useBackproj)
        {
            //SetStartTime();
            stepSize = LSBackproject(vertGradients, stepSize, lu, treeRoot);
            //Debug.Log("Backprojection complete " + DurationString());
        }
        #endregion

        float length2 = TotalLength();
        float energy2 = CurrentEnergy(treeRoot);
        Debug.Log("Length: " + length1 + " -> " + length2 + ", Energy: " + energy1 + " -> " + energy2);

        lastStepSize = stepSize;
        //Debug.Log(soboDot);
        normZero = soboDot < float.NegativeInfinity;//1e-4f; Never true

        Debug.Log("Step duration: " + DurationString());
        return Tuple.Create(stepSize > lsStepThreshold, true);
    }

    /// <summary>
    /// Add all gradients together: From Vectors, Obstacles and Potentials.
    /// Note that the gradients always point in the direction of MORE Energy.
    /// </summary>
    void AddAllGradients(BVHNode2D treeRoot, Matrix<float> vertGradients)
    {
        if (Repulsion)
        {
            if (treeRoot != null)
                TPE.GetInstance.FillGradientVectorBH(this, treeRoot, vertGradients);
            else
                TPE.GetInstance.FillGradientVectorDirect(this, vertGradients);
        }

        // Add gradient contributions from obstacles
        foreach (Obstacle obs in obstacles)
        {
            if (obs.IsEnabled)
                obs.AddGradient(this, vertGradients);
        }

        //Debug.Log(vertGradients);

        // Add gradient contributions from potentials
        foreach (Potential pot in potentials)
        {
            if (!(pot.delayed && scalingLength))
                pot.AddGradient(vertGradients);
        }

        //Debug.Log(vertGradients);

        float x = 0, y = 0;
        for (int i = 0; i < vertGradients.RowCount; i++)
        {
            x += vertGradients.Row(i)[0];
            y += vertGradients.Row(i)[1];
        }
        //Debug.Log("For the gradient to be correct, the following values must be 0:");
        //Debug.Log(x);
        //Debug.Log(y);
    }

    /// <summary>
    /// A line search algorithm. Moves all vertex positions in the direction of their gradient until the Armijo condition is satisfied.
    /// Afterwards, if the total step size is too small, the original positions are restored.
    /// </summary>
    float LineSearchStep(Matrix<float> gradient, float gradDot, BVHNode2D treeRoot, bool resetStep = false)
    {
        float gradNorm = (float)gradient.FrobeniusNorm();
        float initGuess = (gradNorm > 1) ? 1.0f / gradNorm : 1.0f / Mathf.Sqrt(gradNorm);

        // Use the step size from the previous iteration, if it exists
        if (!resetStep && lastStepSize > Mathf.Max(lsStepThreshold, 1e-5f))
            initGuess = Mathf.Min(lastStepSize * 1.5f, initGuess * 4f); // Making sure our initial Guess isn't too big

        float delta = initGuess;

        SaveCurrentPositions();

        float initialEnergy = CurrentEnergy(treeRoot);
        int numBacktracks = 0, numDoubles = 0;
        int doublingLimit = 0;
        float sigma = 0.01f;
        float newEnergy = initialEnergy;

        if (gradNorm < 1e-10)
            return 0;

        while (delta > lsStepThreshold)
        {
            SetGradientStep(gradient, delta);
            newEnergy = CurrentEnergy(treeRoot);

            float decrease = initialEnergy - newEnergy; // How much has the energy decreased?
            float targetDecrease = sigma * delta * gradNorm * gradDot; // How much do we want to decrease? Armijo-Condition

            // See MathematicalConcepts > WolfeConditions.txt
            // If the energy hasn't been decreased enough to meet the Armijo condition, halve the step size
            if (decrease < targetDecrease)
            {
                delta /= 2f;
                numBacktracks++;
            }
            else if (decrease >= targetDecrease && numBacktracks == 0 && numDoubles < doublingLimit)
            {
                delta *= 2f;
                numDoubles++;
            }
            // Otherwise, accept the current step
            else
            {
                SetGradientStep(gradient, delta);
                break;
            }
        }

        if (delta <= lsStepThreshold)   // Restore initial positions if step size goes to 0
        {
            Debug.Log("delta less than threshold");
            RestoreOriginalPositions();
            return 0;
        }

        return delta;
    }

    /// <summary>
    /// Computes the Sobolev-Slobodeckij Gradient. (In rep-curves, this method's name is "ProjectGradient")
    /// </summary>
    float ProjectSoboSloboGradient(Matrix<float> gradients, ref LU<float> sps1_lu)
    {
        Matrix<float> l2gradients = Matrix<float>.Build.DenseOfMatrix(gradients);
        Matrix<float> sps1 = null; // left matrix of saddle point system
        //Debug.Log(l2gradients);

        // Assemble the Sobolev gram matrix with constraints
        SoboSlobo.GetInstance.Sobolev3XWithConstraints(this, ref sps1);

        // Factorize and solve
        sps1_lu = sps1.LU();
        ProjectGradient(sps1_lu, gradients);

        float soboDot = 0;

        //Debug.Log(gradients);
        //Debug.Log(l2gradients);

        for (int i = 0; i < gradients.RowCount; i++)
        {
            soboDot += Vector2.Dot(CurveGenUtils.SelectRow(l2gradients, i), CurveGenUtils.SelectRow(gradients, i));
        }

        return soboDot;
    }

    /// <summary>
    /// (In rep-curves, this method's name is "ProjectSoboSloboGradient")
    /// </summary>
    float ProjectGradient(LU<float> sps1_lu, Matrix<float> gradients)
    {
        // If using per-edge length constraints, then the matrix has all coordinates merged, so we only need one solve
        int numVerts = NumVerts();
        Vector<float> sps3 = Vector<float>.Build.Dense(constraintSet.NumConstraintRows() + constraintSet.NumExpectedCols());

        // Fill in RHS (right hand side) with all coordinates
        for (int i = 0; i < numVerts; i++)
        {
            sps3[2 * i] = gradients[i, 0];
            sps3[2 * i + 1] = gradients[i, 1];
        }

        // Solve for all coordinates
        Vector<float> sps2 = sps1_lu.Solve(sps3);

        for (int i = 0; i < numVerts; i++)
        {
            CurveGenUtils.SetRow(gradients, i, new Vector2(sps2[2 * i], sps2[2 * i + 1]));
        }

        return 1;
    }

    /// <summary>
    /// Saves the current vertex positions, in case we need to revert
    /// </summary>
    void SaveCurrentPositions()
    {
        originalPositions = positions.Clone();
    }

    void RestoreOriginalPositions()
    {
        positions = originalPositions;
    }

    /// <summary>
    /// Moves all vertices by their gradient times delta
    /// </summary>
    void SetGradientStep(Matrix<float> gradient, float delta)
    {
        positions = originalPositions - delta * gradient;
    }

    /// <summary>
    /// Returns the current energy of the whole curve
    /// </summary>
    public float CurrentEnergy(BVHNode2D treeRoot)
    {
        float energy = 0;

        if (Repulsion)
        {
            if (treeRoot != null) energy += EnergyBH(treeRoot);
            else energy += TPE.GetInstance.TpeTotal(this);
        }

        foreach (Obstacle obs in obstacles)
        {
            if (obs.IsEnabled)
                energy += obs.ComputeEnergy(this);
        }

        foreach (Potential pot in potentials)
        {
            if (!(pot.delayed && scalingLength))
                energy += pot.CurrentValue();
        }

        return energy;
    }

    private float EnergyBH(BVHNode2D treeRoot)
    {
        return treeRoot.TotalEnergy(this, treeRoot);
    }

    #endregion

    #region Constraints and Backprojection

    private float LSBackproject(Matrix<float> gradients, float initGuess, LU<float> lu, BVHNode2D root)
    {
        float delta = initGuess;
        int attempts = 0;

        while ((delta > lsStepThreshold || UsingConstraint(ConstraintType.Length)) && attempts < 10)
        {
            attempts++;
            SetGradientStep(gradients, delta);
            if (root != null)
                // Update the centers of mass to reflect the new positions
                root.RecomputeCentersOfMass(this);

            for (int i = 0; i < 3; i++)
            {
                float maxValue = BackprojectConstraints(lu);
                if (maxValue < backproj_threshold)
                {
                    //Debug.Log("Backprojection successful after " + attempts + " attempts");
                    //Debug.Log("Used " + (i + 1) + " Newton steps on successful attempt");
                    return delta;
                }
            }

            delta *= 0.5f;
        }
        Debug.Log("Couldn't make backprojection succeed after " + attempts + " attempts (initial step " + initGuess + ")");
        SetGradientStep(gradients, 0);
        BackprojectConstraints(lu);
        return delta;
    }

    private float BackprojectConstraints(LU<float> lu)
    {
        int numVerts = NumVerts();
        Vector<float> b = Vector<float>.Build.Dense(constraintSet.NumExpectedCols() + constraintSet.NumConstraintRows());

        // If using per-edge length constraints, matrix has all coordinates merged, so we only need one solve

        // Fill RHS with negative constraint values
        float maxViolation = constraintSet.FillConstraintValues(b, constraintTargets, 2 * numVerts);

        // Solve for correction
        Vector<float> corr = lu.Solve(b);

        // Apply correction
        for (int i = 0; i < numVerts; i++)
        {
            Vector2 correction = new Vector2(corr[2 * i], corr[2 * i + 1]);
            CurveVertex pt = verts[i];
            pt.SetPosition(pt.Position() + correction);
            //if (i == 0)
            //    Debug.Log(pt.Position());
        }

        // Compute constraint violation after correction
        maxViolation = constraintSet.FillConstraintValues(b, constraintTargets, 2 * numVerts);
        //Debug.Log("Constraint value: " + maxViolation);

        return maxViolation;
    }

    private void MoveLengthTowardsTarget()
    {
        if (!UsingConstraint(ConstraintType.Length))
            return;

        float currLength = TotalLength();
        int cStart = constraintSet.StartIndexOfConstraint(ConstraintType.Length);
        int numRows = constraintSet.NumRowsForConstraint(ConstraintType.Length);

        float diff = targetLength - currLength;
        float toAdd = diff > 0 ? lengthScaleStep : -lengthScaleStep;
        if (Mathf.Abs(diff) < lengthScaleStep * numRows)
            toAdd = diff / numRows;

        //Debug.Log("constraintTargets before: " + constraintTargets);

        for (int i = cStart; i < cStart + numRows; i++)
            constraintTargets[i] += toAdd;

        //Debug.Log("constraintTargets after: " + constraintTargets);
    }

    public bool TargetLengthReached()
    {
        float currLength = TotalLength();
        return Mathf.Abs(currLength - targetLength) <= lengthScaleStep;
    }

    #endregion

    #region Helper Functions

    internal int NumPins()
    {
        throw new NotImplementedException();
    }

    internal static int NumTangentPins()
    {
        throw new NotImplementedException();
    }

    internal int NumPinnedToSurface()
    {
        throw new NotImplementedException();
    }

    public bool UsingConstraint(ConstraintType type)
    {
        return constraintSet.ConstraintTypes.Contains(type);
    }

    internal Vector2 Barycenter()
    {
        Vector2 center = Vector2.zero;
        float totalMass = 0;
        int numVerts = NumVerts();

        for (int i = 0; i < numVerts; i++)
        {
            float m = verts[i].AvgLength();
            center += verts[i].Position() * m;
            totalMass += m;
        }

        return center / totalMass;
    }

    #endregion

    /// <summary>
    /// Subdivides a curve. Only supports a curve where every vertex has exactly two edges
    /// and every edge exactly two vertices. Curve can be open or closed.
    /// </summary>
    /// <returns></returns>
    public EnergyCurve Subdivide()
    {
        Debug.Log("Subdividing");

        int numEdges = NumEdges();

        List<Vector2> newPositions = new List<Vector2>();

        for (int i = 0; i < numEdges; i++)
        {
            CurveEdge e_i = edges[i];

            CurveVertex prev = e_i.GetPrevVertex();
            CurveVertex next = e_i.GetNextVertex();

            if (i == 0) newPositions.Add(prev.Position());

            Vector2 midpoint = (prev.Position() + next.Position()) * 0.5f;
            newPositions.Add(midpoint);

            if (i < numEdges - 1 || !curveClosed) // If open, always add the last vertex. If closed, only add if it isn't the same as the start vertex.
                newPositions.Add(next.Position());
        }

        //AddVertsAtPositions(newPositions);

        InitVertsEdgesFromPositions(newPositions);

        InitConstraints(constraintSet.ConstraintTypes);

        // ToDo: Replace pinned vertices

        return this;
    }

    /// <summary>
    /// Does everything that should be done after sclaing. Returns true iff the repulsion should stop afterwards.
    /// </summary>
    /// <returns></returns>
    bool AfterScaling()
    {
        //foreach (Potential pot in potentials)
        //{
        //    if (pot.delayed == true)
        //    {
        //        Debug.Log("Target Length reached and detected delayed potential: Resetting Line Search Step");
        //        resetStepSize = true;
        //    }
        //}

        if (rotateAfterScaling)
        {
            Debug.Log(positions);

            float getAvgAngle()
            {
                var angles = new List<float>();

                foreach (CurveEdge edge in edges)
                {
                    var dir = edge.Tangent();

                    if (dir.x < 0 && dir.y >= 0) dir = new Vector2(dir.y, -dir.x);  // Quadrant 2
                    if (dir.x < 0 && dir.y < 0) dir = -dir;                         // Quadrant 3
                    if (dir.x >= 0 && dir.y < 0) dir = new Vector2(-dir.y, dir.x);  // Quadrant 4

                    Vector2 refDir = Vector2.up;
                    if (dir.x > dir.y) refDir = Vector2.right;

                    float angle = Vector2.SignedAngle(dir, refDir);
                    //Debug.Log(angle);
                    angle *= Mathf.Deg2Rad;
                    angles.Add(angle);
                }

                return angles.Average();
            }


            Debug.Log(getAvgAngle());

            float rotAngle = -getAvgAngle();

            float[,] rotArray = { { Mathf.Cos(rotAngle), -Mathf.Sin(rotAngle) }, { Mathf.Sin(rotAngle), Mathf.Cos(rotAngle) } };
            Matrix<float> rotMatrix = Matrix<float>.Build.DenseOfArray(rotArray);
            positions *= rotMatrix;

            Debug.Log(rotMatrix);
            Debug.Log(positions);
            Debug.Log(getAvgAngle());

            return true;
        }

        return false;
    }

    public List<Vector2> Polyline => verts.Select(x => x.Position()).ToList();

    #region Serialization

    public void S_TryInitObstacles()
    {
        obstacles = new();
        foreach (var s_obs in s_obstacles)
            obstacles.Add(new(p_exp: s_obs.p_exp, weight: s_obs.weight, numPoints: s_obs.numPoints, radius: s_obs.radius, center: s_obs.center));
    }

    public void S_SerializeObstaclePositions()
    {
        s_obstacles = new();
        foreach (var obs in obstacles)
            s_obstacles.Add(new S_ObstacleConfig { p_exp = obs.p_exp, weight = obs.weight, numPoints = obs.numPoints, radius = obs.radius, center = obs.center });
    }

    public void S_TryInitConstraints()
    {

    }

    public void S_SerializeConstraintTypes()
    {

    }

    #endregion
}

[Serializable]
public enum GenMode
{
    Bezier, Circular
}

[Serializable]
public enum RepulsionType
{
    Normal, Sobolev
}

[Serializable]
public class EnergyCurve_EditorConfig // ToDo: Outsource deacObsAfterScaling, obstacles and numObstacles into separate class
{
    public bool curveClosed;
    public float targetLengthScale;
    public bool deacObsAfterScaling;
    public bool rotateAfterScaling;
    public bool noRepulsionAfterScaling;
    public GenModeConfig genModeConfig;
    public List<ObstacleConfig> obstacles;
    public int numObstacles;
    public List<PotentialConfig> potentials;
    public List<ConstraintType> constraints;

    /// <summary>
    /// 
    /// </summary>
    /// <param name="curveClosed"></param>
    /// <param name="targetLengthScale"></param>
    /// <param name="deacObsAfterScaling"></param>
    /// <param name="genModeConfig"></param>
    /// <param name="obstacles">
    /// If null: Not overriding, should be created algorithmically
    /// If 0 or N elements: Overriding
    /// </param>
    /// <param name="potentials">
    /// Needed for when obstacles should be overriden
    /// </param>
    /// <param name="constraints"></param>
    public EnergyCurve_EditorConfig(bool curveClosed, float targetLengthScale, bool deacObsAfterScaling, bool rotateAfterScaling, bool noRepulsionAfterScaling, GenModeConfig genModeConfig, List<ObstacleConfig> obstacles, int numObstacles, List<PotentialConfig> potentials, List<ConstraintType> constraints)
    {
        this.curveClosed = curveClosed;
        this.targetLengthScale = targetLengthScale;
        this.deacObsAfterScaling = deacObsAfterScaling;
        this.rotateAfterScaling = rotateAfterScaling;
        this.noRepulsionAfterScaling = noRepulsionAfterScaling;
        this.genModeConfig = genModeConfig;
        this.obstacles = obstacles;
        this.numObstacles = numObstacles;
        this.potentials = potentials;
        this.constraints = constraints;
    }
}

[Serializable]
public class ObstacleConfig
{
    public float weight, radius;
    public int numPoints;
    public Vector2 center;

    public ObstacleConfig(float weight, float radius, int numPoints, Vector2 center)
    {
        this.weight = weight;
        this.radius = radius;
        this.numPoints = numPoints;
        this.center = center;
    }
}

[Serializable]
public class PotentialConfig
{
    public PotentialType type;
    public float weight;
    public bool delayed;

    public PotentialConfig(PotentialType type, float weight, bool delayed)
    {
        this.type = type;
        this.weight = weight;
        this.delayed = delayed;
    }
}

[Serializable]
public class GenModeConfig
{

}

[Serializable]
public class GenModeConfig_Bezier : GenModeConfig
{
    public int numPoints = 5;
    public float spacing = 1;
    public float radius = 10;
}

[Serializable]
public class GenModeConfig_Circular : GenModeConfig
{
    public int numPoints = 50;
    public float radius = 10;
}

[Serializable]
struct S_ObstacleConfig
{
    public float p_exp;
    public float weight;
    public int numPoints;
    public float radius;
    public Vector2 center;
}