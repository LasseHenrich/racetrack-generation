using MathNet.Numerics.LinearAlgebra;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshObstacle : Obstacle
{

    public override void AddGradient(EnergyCurve curve, Matrix<float> gradient)
    {
        throw new System.NotImplementedException();
    }

    public override float ComputeEnergy(EnergyCurve curve)
    {
        throw new System.NotImplementedException();
    }
}
