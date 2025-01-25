using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEditor;
using UnityEngine;

public class ToolWindow_PlayMode : MonoBehaviour
{
    GUIStyle style_Label;
    GUIStyle style_InlineLabel;

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
        style_Label = new GUIStyle(GUI.skin.label);
        style_InlineLabel = new GUIStyle(style_Label) { alignment = TextAnchor.MiddleLeft };
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
        GUILayout.BeginHorizontal();
        GUILayout.Label(name, style_InlineLabel);
        input = GUILayout.TextField(input);
        GUILayout.EndHorizontal();
    }

    protected void CreateIntField(string name, ref int input)
    {
        GUILayout.BeginHorizontal();
        GUILayout.Label(name, style_InlineLabel);
        input = int.TryParse(GUILayout.TextField(input.ToString()), out int result) ? result : input;
        GUILayout.EndHorizontal();
    }
    protected void CreateFloatField(string name, ref float input)
    {
        GUILayout.BeginHorizontal();
        GUILayout.Label(name, style_InlineLabel);
        input = float.TryParse(GUILayout.TextField(input.ToString()), out float result) ? result : input;
        GUILayout.EndHorizontal();
    }

    protected void CreateVector2Field(string name, ref Vector2 input)
    {
        CreateFloatField(name + " X", ref input.x);
        CreateFloatField(name + " Y", ref input.y);
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
        GUILayout.Label($"{text}: {value}");
        float rawValue = GUILayout.HorizontalSlider(value, start, end);
        value = Mathf.RoundToInt(rawValue);
    }

    protected void CreateSliderFloat(ref float value, string text, float start, float end, float snap = 0.5f)
    {
        GUILayout.Label($"{text}: {value:F1}");
        float rawValue = GUILayout.HorizontalSlider(value, start, end);
        value = Mathf.Round(rawValue / snap) * snap;
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

    protected void CreateEnumSelection<T>(string name, T currentSelection, Action<T> onValueChanged) where T : Enum
    {
        // Create a section with the enum name
        GUILayout.Label(name, style_Label);

        // Get all values of the enum
        T[] enumValues = (T[])Enum.GetValues(typeof(T));

        // Create a button for each enum value
        foreach (T enumValue in enumValues)
        {
            // Determine if this is the currently selected value
            bool isSelected = currentSelection.Equals(enumValue);

            // Create a button with custom styling based on selection
            GUIStyle buttonStyle = new GUIStyle(GUI.skin.button)
            {
                fontStyle = isSelected ? FontStyle.Bold : FontStyle.Normal,
                normal = {
                textColor = isSelected ? Color.green : Color.white
            }
            };

            // Create the button
            if (GUILayout.Button(enumValue.ToString(), buttonStyle))
            {
                onValueChanged(enumValue);
            }
        }
    }

    #endregion
}
