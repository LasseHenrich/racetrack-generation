using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;

public class CameraController : MonoBehaviour
{
    [SerializeField] float mouseSensitivity = 100f;
    [SerializeField] float moveSpeed = 10f;

    bool isCamMovementActive = false;

    void Update()
    {

        if (Input.GetKeyDown(KeyCode.Space))
        {
            isCamMovementActive = !isCamMovementActive;
        }

        if (!isCamMovementActive)
            return;

        float mouseX = Input.GetAxis("Mouse X") * mouseSensitivity * Time.deltaTime;
        float mouseY = Input.GetAxis("Mouse Y") * mouseSensitivity * Time.deltaTime;

        // Rotate camera based on mouse X and Y movement
        transform.Rotate(Vector3.up * mouseX);
        transform.Rotate(Vector3.left * mouseY);

        // Movement
        float x = Input.GetAxis("Horizontal") * moveSpeed * Time.deltaTime;
        float z = Input.GetAxis("Vertical") * moveSpeed * Time.deltaTime;

        transform.Translate(x, 0, z);
    }
}
