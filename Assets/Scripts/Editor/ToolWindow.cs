using System;
using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class ToolWindow : EditorWindow
{
    protected virtual string ToolLabel { get; }

    Vector2 scrollPos;

    GUIStyle style_Section;
    GUIStyle style_Label;

    Dictionary<int, bool> foldouts;

    // Start is called before the first frame update
    protected virtual void OnEnable()
    {
        foldouts = new();
    }

    protected virtual void OnGUI()
    {
        style_Section = new GUIStyle(EditorStyles.foldout)
        {
            fontStyle = FontStyle.Bold,
            alignment = TextAnchor.MiddleLeft,
            padding = new RectOffset(15, 0, 0, 0),
        };

        style_Label = new GUIStyle(EditorStyles.label);

        GUILayout.BeginVertical();
        GUILayout.Space(10);
        scrollPos = EditorGUILayout.BeginScrollView(scrollPos);
    }

    #region GUI Methods

    protected void CreateSection(string name, Action content)
    {
        CreateFoldout(name, content, 20, 6);
    }
    protected void CreateSubSection(string name, Action content)
    {
        CreateFoldout(name, content, 16, 4);
    }

    protected void CreateSubSubSection(string name, Action content)
    {
        CreateFoldout(name, content, 14, 3);
    }

    protected void CreateLabel(string name, int fontSize = 12)
    {
        GUIStyle style = new(style_Label) { fontSize = fontSize };
        GUILayout.Label(name, style);
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

    protected void CreateFoldout(string name, Action content, int fontSize, int bottomSpacing)
    {
        int foldoutHash = content.GetHashCode();

        GUIStyle style = new(style_Section) { fontSize = fontSize };

        bool initialization = false;
        if (!foldouts.ContainsKey(foldoutHash))
        {
            initialization = true;
            foldouts.Add(foldoutHash, true); // Initially everything must be folded out, as otherwise some initializations fail
        }

        foldouts[foldoutHash] = EditorGUILayout.Foldout(
            foldout: foldouts.ContainsKey(foldoutHash) && foldouts[foldoutHash],
            content: name,
            toggleOnLabelClick: true,
            style: style
        );

        if (foldouts[foldoutHash])
        {
            GUILayout.BeginHorizontal();
            GUILayout.Space(10);

            GUILayout.BeginVertical();
            GUILayout.Space(bottomSpacing);

            content();

            GUILayout.Space(bottomSpacing); // 2. Additional spacing under last element
            GUILayout.EndVertical();
            GUILayout.EndHorizontal();
        }
        GUILayout.Space(bottomSpacing); // 1. Default bottom spacing when collapsed

        if (initialization)
        {
            foldouts[foldoutHash] = false;
        }
    }

    #endregion
}
