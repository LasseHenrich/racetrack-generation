using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public abstract class Obstacle
{
    bool enabled;
    public bool IsEnabled { get { return enabled; } }

    public Obstacle()
    {
        enabled = true;
    }

    public void Disable()
    {
        enabled = false;
    }

    public void Enable()
    {
        enabled = true;
    }

    public abstract void AddGradient(EnergyCurve curve, Matrix<float> gradient);
    public abstract float ComputeEnergy(EnergyCurve curve);
}
