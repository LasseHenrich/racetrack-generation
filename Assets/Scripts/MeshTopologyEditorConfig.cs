using System.Collections.Generic;
using System;
using UnityEngine;

[Serializable]
public class MeshTopologyEditorConfig
{
    public Mesh mesh;
    public int numMaterials = 3;
    public List<Material> materials = new();
    public bool showFoldout;

    public MeshTopologyEditorConfig(Mesh mesh, int numMaterials, List<Material> materials, bool showFoldout)
    {
        this.mesh = mesh;
        this.numMaterials = numMaterials;
        this.materials = materials;
        this.showFoldout = showFoldout;
    }
}