using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class ToolWindow : EditorWindow
{
    protected virtual string ToolLabel { get; }

    Vector2 scrollPos;

    GUIStyle style_Section;
    GUIStyle style_SubSection;
    GUIStyle style_SubSubSection;
    GUIStyle style_Label;

    // Start is called before the first frame update
    protected virtual void OnEnable()
    {

        Debug.Log("OnEnable");
    }

    protected virtual void OnGUI()
    {
        style_Section = new GUIStyle(EditorStyles.boldLabel)
        {
            fontSize = 24,
            fixedHeight = 42
        };

        style_SubSection = new GUIStyle(EditorStyles.boldLabel)
        {
            fontSize = 20,
            fixedHeight = 36
        };

        style_SubSubSection = new GUIStyle(EditorStyles.boldLabel)
        {
            fontSize = 16,
            fixedHeight = 30
        };

        style_Label = new GUIStyle(EditorStyles.label);

        GUILayout.BeginVertical();
        scrollPos = EditorGUILayout.BeginScrollView(scrollPos);

        CreateSection(ToolLabel);
    }

    #region GUI Methods

    protected void CreateSection(string name)
    {
      
        GUILayout.Label(name, style_Section);
    }
    protected void CreateSubSection(string name)
    {
        GUILayout.Label(name, style_SubSection);
    }

    protected void CreateSubSubSection(string name)
    {
        GUILayout.Label(name, style_SubSubSection);
    }

    protected void CreateLabel(string name)
    {
        GUILayout.Label(name, style_Label);
    }

    protected void CreateTextField(string name, ref string input)
    {
        input = EditorGUILayout.TextField(name, input);
    }

    protected void CreateIntField(string name, ref int input)
    {
     
        input = EditorGUILayout.IntField(name, input);
    }
    protected void CreateFloatField(string name, ref float input)
    {
        input = EditorGUILayout.FloatField(name, input);
    }

    protected void CreateVector2Field(string name, ref Vector2 input)
    {
        input = EditorGUILayout.Vector2Field(name, input);
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
        Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
        position.height *= 0.5f;
        value = (int)EditorGUI.Slider(position, text, value, start, end);

        position.y += position.height;
        position.x += EditorGUIUtility.labelWidth;
        position.width -= EditorGUIUtility.labelWidth + 54;

        GUIStyle style = GUI.skin.label;
        style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
        style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
    }

    protected void CreateSliderFloat(ref float value, string text, float start, float end)
    {
        Rect position = EditorGUILayout.GetControlRect(false, 2 * EditorGUIUtility.singleLineHeight);
        position.height *= 0.5f;
        value = EditorGUI.Slider(position, text, value, start, end);

        position.y += position.height;
        position.x += EditorGUIUtility.labelWidth;
        position.width -= EditorGUIUtility.labelWidth + 54;

        GUIStyle style = GUI.skin.label;
        style.alignment = TextAnchor.UpperLeft; EditorGUI.LabelField(position, start.ToString(), style);
        style.alignment = TextAnchor.UpperRight; EditorGUI.LabelField(position, end.ToString(), style);
    }

    protected void CreateCurveField(string name, ref AnimationCurve curve)
    {
        curve = EditorGUILayout.CurveField(name, curve);
    }

    protected void CreateFoldout(string name, ref bool foldout)
    {
        Rect position = EditorGUILayout.GetControlRect(false, EditorGUIUtility.singleLineHeight);
        foldout = EditorGUI.Foldout(position, foldout, name);
    }

    #endregion
}
