import sys
import os
import argparse
import subprocess
from PyQt5.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QComboBox,
    QMessageBox,
    QCheckBox,
    QProgressDialog,
    QFormLayout,
)
from PyQt5.QtCore import Qt, QTimer, QRegExp
from PyQt5.QtGui import QIntValidator, QRegExpValidator
import re


class BarcodeTargetSeekerGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.parser = self.create_parser()
        self.initUI()

    def create_parser(self):
        parser = argparse.ArgumentParser(
            description="Map barcodes to a circular genome"
        )
        parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
        parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
        parser.add_argument("pam", help="PAM sequence", type=str)
        parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)
        parser.add_argument(
            "--pam_direction",
            choices=["upstream", "downstream"],
            default="downstream",
            help="Direction of the PAM sequence",
        )
        parser.add_argument(
            "--json",
            action="store_true",
            default=False,
            help="Output results in JSON format",
        )
        return parser

    def initUI(self):
        main_layout = QVBoxLayout()

        # Add a back button
        back_btn = QPushButton("← Back to Main Menu")
        back_btn.clicked.connect(
            lambda: self.parent().setCurrentWidget(
                self.parent().widget(0)  # Assuming welcome page is at index 0
            )
        )
        main_layout.addWidget(back_btn)

        # Title
        title = QLabel("Barcode Target Seeker")
        title.setStyleSheet("font-size: 20px; font-weight: bold; margin: 10px 0;")
        main_layout.addWidget(title)

        # Dynamic argument input fields
        self.argument_widgets = {}
        self.optional_widgets = {}

        # Separate positional and optional arguments
        positional_actions = [
            action
            for action in self.parser._actions
            if action.dest not in ["help", None, "pam_direction", "json"]
        ]
        optional_actions = [
            action
            for action in self.parser._actions
            if action.dest in ["pam_direction", "json"]
        ]

        # Use QFormLayout for inputs
        form_layout = QFormLayout()

        # Positional arguments
        for action in positional_actions:
            if action.dest == "sgrna_file":
                label_text = "sgRNA File"
            else:
                label_text = action.dest.replace("_", " ").title()
            label = QLabel(label_text)
            label.setMinimumWidth(120)

            if action.dest in ["sgrna_file", "genome_file"]:
                input_widget = QLineEdit()
                browse_button = QPushButton("Browse")
                browse_button.clicked.connect(
                    lambda checked, dest=action.dest: self.browse_file(dest)
                )
                file_layout = QHBoxLayout()
                file_layout.addWidget(input_widget)
                file_layout.addWidget(browse_button)
                form_layout.addRow(label, file_layout)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "mismatches":
                input_widget = QComboBox()
                input_widget.addItems(["0", "1", "2"])
                input_widget.setCurrentText("0")  # Set default to 0
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

            elif action.dest == "pam":
                input_widget = QLineEdit("NGG")  # Set default to NGG
                # Add validator for PAM sequence (only ACGTN letters)
                validator = QRegExpValidator(QRegExp("^[ACGTNacgtn]+$"))
                input_widget.setValidator(validator)
                form_layout.addRow(label, input_widget)
                self.argument_widgets[action.dest] = input_widget

        main_layout.addLayout(form_layout)

        # Optional arguments
        optional_label = QLabel("Optional Arguments")
        optional_label.setStyleSheet("font-weight: bold; margin-top: 15px;")
        main_layout.addWidget(optional_label)

        for action in optional_actions:
            arg_layout = QHBoxLayout()
            label = QLabel(action.dest.replace("_", " ").title())
            label.setMinimumWidth(120)

            if action.dest == "pam_direction":
                input_widget = QComboBox()
                input_widget.addItems(action.choices)
                input_widget.setCurrentText("downstream")  # default
                arg_layout.addWidget(label)
                arg_layout.addWidget(input_widget)
                self.optional_widgets[action.dest] = input_widget

            elif action.dest == "json":
                input_widget = QCheckBox("Output as JSON")
                arg_layout.addWidget(input_widget)
                self.optional_widgets[action.dest] = input_widget

            main_layout.addLayout(arg_layout)

        # Add input field for output filename
        self.output_filename_label = QLabel("Output Filename:")
        self.output_filename_input = QLineEdit()
        self.output_filename_input.setPlaceholderText("Enter output filename")
        main_layout.addWidget(self.output_filename_label)
        main_layout.addWidget(self.output_filename_input)

        # Run button
        run_button = QPushButton("Run Barcode Target Seeker")
        run_button.setStyleSheet(
            """
            QPushButton {
                background-color: #007bff;
                color: white;
                padding: 10px;
                border-radius: 4px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #0056b3;
            }
        """
        )
        run_button.clicked.connect(self.run_script)
        main_layout.addWidget(run_button)

        # Add stretch to push everything to the top
        main_layout.addStretch()

        self.setLayout(main_layout)

    def browse_file(self, dest):
        if dest == "sgrna_file":
            filter_str = "FASTA Files (*.fasta *.fa);;All Files (*)"
        elif dest == "genome_file":
            filter_str = "GenBank Files (*.gb *.gbk);;All Files (*)"
        else:
            filter_str = "All Files (*)"
        filename, _ = QFileDialog.getOpenFileName(
            self,
            f"Select {self.argument_widgets[dest].placeholderText()}",
            "",
            filter_str,
        )
        if filename:
            self.argument_widgets[dest].setText(filename)

    def run_script(self):
        # Validate required inputs
        args = []
        for dest, widget in self.argument_widgets.items():
            if isinstance(widget, QComboBox):
                args.append(widget.currentText())
            elif not widget.text():
                QMessageBox.warning(
                    self, "Input Error", f"{dest.replace('_', ' ').title()} is required"
                )
                return
            else:
                args.append(widget.text())

        # Add optional arguments
        for dest, widget in self.optional_widgets.items():
            if dest == "pam_direction":
                if widget.currentText() != "downstream":
                    args.extend(["--pam_direction", widget.currentText()])
            elif dest == "json":
                if widget.isChecked():
                    args.append("--json")

        # Get output filename and enforce it
        output_filename = self.output_filename_input.text().strip()
        if not output_filename:
            QMessageBox.warning(self, "Input Error", "Output Filename is required.")
            return

        # Enforce file extension based on JSON checkbox
        if self.optional_widgets["json"].isChecked():
            if not output_filename.lower().endswith(".json"):
                output_filename += ".json"
        else:
            if not output_filename.lower().endswith(".tsv"):
                output_filename += ".tsv"

        # Add validation for 'pam' input
        pam_input = self.argument_widgets["pam"].text().upper().strip()
        if not pam_input:
            QMessageBox.warning(self, "Input Error", "PAM sequence is required.")
            return
        elif not re.match("^[ACGTN]+$", pam_input):
            QMessageBox.warning(
                self, "Input Error", "Invalid PAM sequence. Use A, C, G, T, or N."
            )
            return

        # Create progress dialog
        progress = QProgressDialog(
            "Running Barcode Target Seeker...", "Cancel", 0, 0, self
        )
        progress.setWindowModality(Qt.WindowModal)
        progress.setAutoReset(False)
        progress.setAutoClose(False)
        progress.show()

        # Construct the full command to run the script
        script_path = os.path.join(os.path.dirname(__file__), "targets.py")
        full_command = ["python", script_path] + args

        try:
            # Run the script as a subprocess and redirect stdout to the file
            with open(output_filename, "w") as f_out:
                process = subprocess.Popen(
                    full_command,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                )

            # Create a timer to check process status
            timer = QTimer(self)
            timer.timeout.connect(lambda: self.check_process(process, progress, timer))
            timer.start(100)  # Check every 100ms

        except Exception as e:
            progress.close()
            QMessageBox.critical(self, "Error", str(e))

    def check_process(self, process, progress, timer):
        if process.poll() is not None:
            # Process finished
            timer.stop()
            progress.close()

            stderr = process.stderr.read()

            if process.returncode == 0:
                # Success: show output in a message box
                output_dialog = QMessageBox(self)
                output_dialog.setWindowTitle("Script Output")
                output_dialog.setText("Script executed successfully!")
                output_dialog.setStandardButtons(QMessageBox.Ok)
                output_dialog.exec_()
            else:
                # Error: show error message
                error_dialog = QMessageBox(self)
                error_dialog.setIcon(QMessageBox.Critical)
                error_dialog.setWindowTitle("Error")
                error_dialog.setText("Script execution failed")
                error_dialog.setDetailedText(stderr)
                error_dialog.setStandardButtons(QMessageBox.Ok)
                error_dialog.exec_()
