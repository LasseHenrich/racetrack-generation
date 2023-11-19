using System;
using System.Collections.Generic;
using UnityEngine;

public class ToolWindow : MonoBehaviour
{
    GUIStyle style_Label;

    Dictionary<int, bool> foldouts;

    // Start is called before the first frame update
    protected virtual void OnEnable()
    {
        useGUILayout = true;
        Debug.Log("OnEnable");
        foldouts = new();
    }

    protected virtual void OnGUI()
    {
        style_Label = new GUIStyle();
    }

    #region GUI Methods

    protected void CreateSection(string name, Action content)
    {
        CreateFoldout(name, content, 24);
    }
    protected void CreateSubSection(string name, Action content)
    {
        CreateFoldout(name, content, 20);
    }

    protected void CreateSubSubSection(string name, Action content)
    {
        CreateFoldout(name, content, 16);
    }

    protected void CreateLabel(string name)
    {
        GUILayout.Label(name, style_Label);
    }

    protected void CreateTextField(string name, ref string input)
    {
        input = GUILayout.TextField(name, input);
    }

    protected void CreateIntField(string name, ref int input)
    {
        input = int.TryParse(GUILayout.TextField(name + ": " + input), out int result) ? result : input;
    }
    protected void CreateFloatField(string name, ref float input)
    {
        input = float.TryParse(GUILayout.TextField(name + ": " + input), out float result) ? result : input;
    }

    protected void CreateVector2Field(string name, ref Vector2 input)
    {
        CreateFloatField(name + " X", ref input.x);
        CreateFloatField(name + " Y", ref input.y);
    }

    protected void CreateVector3Field(string name, ref Vector3 input)
    {
        CreateFloatField(name + " X", ref input.x);
        CreateFloatField(name + " Y", ref input.y);
        CreateFloatField(name + " Z", ref input.z);
    }

    protected void CreateButton(string name, System.Action callback)
    {
        if (GUILayout.Button(name)) callback();
    }

    protected void CreateCheckbox(string name, ref bool trigger)
    {
        trigger = GUILayout.Toggle(trigger, name);
    }

    protected void CreateCheckbox_Dict<T>(string name, Dictionary<T, bool> dict, T key)
    {
        dict[key] = GUILayout.Toggle(dict[key], name);
    }

    protected void CreateSlider(ref int value, string text, int start, int end)
    {
        /*
        Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
        position.height *= 0.5f;
        value = (int)EditorGUI.Slider(position, text, value, start, end);

        position.y += position.height;
        position.x += EditorGUIUtility.labelWidth;
        position.width -= EditorGUIUtility.labelWidth + 54;

        GUIStyle style = GUI.skin.label;
        style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
        style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
        */
        value = (int)GUILayout.HorizontalSlider(value, start, end);
    }

    protected void CreateSliderFloat(ref float value, string text, float start, float end)
    {
        /*
        Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
        position.height *= 0.5f;
        value = EditorGUI.Slider(position, text, value, start, end);

        position.y += position.height;
        position.x += EditorGUIUtility.labelWidth;
        position.width -= EditorGUIUtility.labelWidth + 54;

        GUIStyle style = GUI.skin.label;
        style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
        style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
        */
        value = GUILayout.HorizontalSlider(value, start, end);
    }

    protected void CreateCurveField(string name, ref AnimationCurve curve)
    {
        throw new NotImplementedException();
        //curve = EditorGUILayout.CurveField(name, curve);
    }

    protected void _CreateFoldout(string name, ref bool foldout)
    {
        /*
        Rect position = EditorGUILayout.GetControlRect(false, EditorGUIUtility.singleLineHeight);
        foldout = EditorGUI.Foldout(position, foldout, name);
        */
        foldout = GUILayout.Toggle(foldout, name);
    }

    protected void CreateFoldout(string name, Action content, int fontSize)
    {
        GUIStyle foldoutStyle = new()
        {
            fontSize = fontSize,
            fontStyle = FontStyle.Bold,
            normal = {
                textColor = Color.white              // Normal text color, for example black (change as needed)
            },
            hover = {
                background = Texture2D.whiteTexture, // A plain white texture
                textColor = Color.blue               // Highlight text color, for example blue (change as needed)
            }
        };


        int foldoutHash = content.GetHashCode();

        /*
        Rect position = EditorGUILayout.GetControlRect(true, EditorGUIUtility.singleLineHeight * 2);
        foldouts[foldoutHash] = EditorGUI.Foldout(
            position,
            foldouts.ContainsKey(foldoutHash) && foldouts[foldoutHash],
            name,
            true,
            foldoutStyle
        );
        */

        foldouts[foldoutHash] = GUILayout.Toggle(
            foldouts.ContainsKey(foldoutHash) && foldouts[foldoutHash],
            name,
            foldoutStyle
        );

        if (foldouts[foldoutHash])
        {
            GUILayout.BeginHorizontal();
            GUILayout.Space(10);
            GUILayout.BeginVertical();

            content();

            GUILayout.EndVertical();
            GUILayout.EndHorizontal();
        }
    }

    protected void CreateEnumSelection<T>(string name, ref T op) where T : Enum
    {
        GUILayout.Label(name);
        string[] enumNames = Enum.GetNames(typeof(T));
        op = (T)(object)GUILayout.SelectionGrid((int)(object)op, enumNames, 2);
    }

    #endregion
}
