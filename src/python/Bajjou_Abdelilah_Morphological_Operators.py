import cv2
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk

class MorphologicalApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Morphological Operations")

        self.static_image = None  # Original image (never changes)
        self.treated_image = None  # Treated image (changes based on operations)
        self.kernel = np.ones((5, 5), np.uint8)
        self.apply_to_treated = tk.BooleanVar(value=False)  # Checkbox state

        # Create widgets
        self.create_widgets()

        # Load placeholder image
        self.load_placeholder()

    def create_widgets(self):
        # Frame for image display
        image_frame = tk.Frame(self.root)
        image_frame.pack(side=tk.TOP, pady=10)

        # Static Image Label
        self.static_label = tk.Label(image_frame, text="Original Image", font=("Arial", 12))
        self.static_label.pack(side=tk.LEFT, padx=10)

        # Treated Image Label
        self.treated_label = tk.Label(image_frame, text="Treated Image", font=("Arial", 12))
        self.treated_label.pack(side=tk.LEFT, padx=10)

        # Frame for controls
        control_frame = tk.Frame(self.root)
        control_frame.pack(side=tk.TOP, pady=10)

        # Load Image Button
        tk.Button(
            control_frame,
            text="Load Image",
            command=self.load_image,
            font=("Arial", 10)
        ).pack(side=tk.LEFT, padx=5)

        # Save Image Button
        tk.Button(
            control_frame,
            text="Save Treated Image",
            command=self.save_image,
            font=("Arial", 10)
        ).pack(side=tk.LEFT, padx=5)

        # 📌 **Reset Image Button**
        tk.Button(
            control_frame,
            text="Reset Treated Image",
            command=self.reset_image,
            font=("Arial", 10)
        ).pack(side=tk.LEFT, padx=5)

        # Checkbox to apply changes to treated image
        self.apply_checkbox = tk.Checkbutton(
            control_frame,
            text="Apply Changes to Treated Image",
            variable=self.apply_to_treated,
            font=("Arial", 10)
        )
        self.apply_checkbox.pack(side=tk.LEFT, padx=5)

        # Threshold Slider
        self.threshold_label = tk.Label(control_frame, text="Threshold: 127", font=("Arial", 10))
        self.threshold_label.pack(side=tk.LEFT, padx=5)

        self.threshold_slider = tk.Scale(
            control_frame,
            from_=0,
            to=255,
            orient=tk.HORIZONTAL,
            command=self.update_threshold,
            font=("Arial", 10)
        )
        self.threshold_slider.set(127)
        self.threshold_slider.pack(side=tk.LEFT, padx=5)

        # Frame for morphological operations buttons
        operations_frame = tk.Frame(self.root)
        operations_frame.pack(side=tk.TOP, pady=10)

        # Morphological Operations Buttons
        operations = [
            "Erode", "Dilate", "Open", "Close",
            "Top Hat (Bright)", "Top Hat (Dark)",
            "Gradient", "Median",
            "Contour (Internal)", "Contour (External)",
            "Morphological Contour"
        ]

        # Organize buttons in a grid (3 rows, 4 columns)
        row, col = 0, 0
        for operation in operations:
            tk.Button(
                operations_frame,
                text=operation,
                command=lambda op=operation: self.apply_operation(op),
                font=("Arial", 10),
                width=20
            ).grid(row=row, column=col, padx=5, pady=5)
            col += 1
            if col > 3:  # 4 buttons per row
                col = 0
                row += 1

    def load_placeholder(self):
        # Create a placeholder image (black with white text)
        placeholder = np.zeros((200, 300), dtype=np.uint8)
        cv2.putText(
            placeholder,
            "No Image Loaded",
            (30, 100),
            cv2.FONT_HERSHEY_SIMPLEX,
            1,
            255,
            2
        )
        self.static_image = placeholder
        self.treated_image = placeholder.copy()
        self.show_image(self.static_image, self.static_label)
        self.show_image(self.treated_image, self.treated_label)

    def load_image(self):
        file_path = filedialog.askopenfilename()
        if file_path:
            self.static_image = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
            self.treated_image = self.static_image.copy()
            self.show_image(self.static_image, self.static_label)
            self.show_image(self.treated_image, self.treated_label)

    def update_threshold(self, value):
        if self.static_image is not None:
            _, thresholded = cv2.threshold(
                self.static_image,
                int(value),
                255,
                cv2.THRESH_BINARY
            )
            self.threshold_label.config(text=f"Threshold: {value}")
            if self.apply_to_treated.get():
                self.treated_image = thresholded
                self.show_image(self.treated_image, self.treated_label)
            else:
                self.show_image(thresholded, self.treated_label)

    def apply_operation(self, operation):
        if self.static_image is not None:
            # Decide which image to apply operations on
            if not self.apply_to_treated.get():
                target_image = self.static_image
            else:
                target_image = self.treated_image

            if operation == "Erode":
                result = cv2.erode(target_image, self.kernel, iterations=1)
            elif operation == "Dilate":
                result = cv2.dilate(target_image, self.kernel, iterations=1)
            elif operation == "Open":
                result = cv2.morphologyEx(target_image, cv2.MORPH_OPEN, self.kernel)
            elif operation == "Close":
                result = cv2.morphologyEx(target_image, cv2.MORPH_CLOSE, self.kernel)
            elif operation == "Top Hat (Bright)":
                result = cv2.morphologyEx(target_image, cv2.MORPH_TOPHAT, self.kernel)
            elif operation == "Top Hat (Dark)":
                result = cv2.morphologyEx(target_image, cv2.MORPH_BLACKHAT, self.kernel)
            elif operation == "Gradient":
                result = cv2.morphologyEx(target_image, cv2.MORPH_GRADIENT, self.kernel)
            elif operation == "Median":
                result = cv2.medianBlur(target_image, 5)
            elif operation == "Contour (Internal)":
                result = target_image - cv2.erode(target_image, self.kernel)
            elif operation == "Contour (External)":
                result = cv2.dilate(target_image, self.kernel) - target_image
            elif operation == "Morphological Contour":
                opening = cv2.morphologyEx(target_image, cv2.MORPH_OPEN, self.kernel)
                closing = cv2.morphologyEx(target_image, cv2.MORPH_CLOSE, self.kernel)
                result = closing - opening
            else:
                result = target_image

            # Update treated image
            self.treated_image = result
            self.show_image(self.treated_image, self.treated_label)

    def save_image(self):
        if self.treated_image is not None:
            file_path = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg *.jpeg")]
            )
            if file_path:
                cv2.imwrite(file_path, self.treated_image)
                messagebox.showinfo("Success", "Treated image saved successfully!")

    # 📌 **Reset Treated Image Method**
    def reset_image(self):
        if self.static_image is not None:
            self.treated_image = self.static_image.copy()
            self.show_image(self.treated_image, self.treated_label)
            messagebox.showinfo("Reset", "Treated image has been reset to the original!")

    def show_image(self, image, label):
        image = Image.fromarray(image)
        image = ImageTk.PhotoImage(image)
        label.config(image=image)
        label.image = image

if __name__ == "__main__":
    root = tk.Tk()
    app = MorphologicalApp(root)
    root.mainloop()
