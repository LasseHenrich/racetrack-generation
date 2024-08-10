using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

[Serializable]
public class SerializedMapObject
{
    [SerializeField]
    public RoadSpline roadSpline;
    [SerializeField]
    public List<TopologyRange> topologies;
    [SerializeField]
    public List<float> fourPercCrossings_1D;  // List<List<float>> not serializable with JSONUtility

    private const string localPath = "Resources/Maps/Generated/SplineAndTops";

    public SerializedMapObject(RoadSpline roadSpline, List<TopologyRange> topologies, List<List<float>> fourPercCrossings)
    {
        this.roadSpline = roadSpline;
        this.topologies = topologies;

        this.fourPercCrossings_1D = new();
        foreach (List<float> crossing in fourPercCrossings)
        {
            this.fourPercCrossings_1D.AddRange(crossing);
        }
    }

    public void SaveToFile(string name)
    {
        string json = JsonUtility.ToJson(this);
        //Debug.Log(json);

        string path = $"{Application.dataPath}/{localPath}/{name}" +
            //$"_{DateTime.Now.Year}{DateTime.Now.Month}{DateTime.Now.Day}{DateTime.Now.Hour}{DateTime.Now.Minute}{DateTime.Now.Second}" +
            $".json";

        using FileStream fs = new(path, FileMode.CreateNew);
        using StreamWriter writer = new(fs);
        writer.Write(json);
    }

    public static SerializedMapObject LoadFromFile(string name)
    {
        string path = $"{Application.dataPath}/{localPath}/{name}.json";

        using FileStream fs = new(path, FileMode.Open);
        using StreamReader reader = new(fs);
        string json = reader.ReadToEnd();

        SerializedMapObject smo = JsonUtility.FromJson<SerializedMapObject>(json);

        smo.roadSpline.controls = new(smo.roadSpline, smo.roadSpline.controls.GetList());

        return smo;
    }
}
