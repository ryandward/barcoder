import sys
import os
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QPushButton,
    QLabel,
    QStackedWidget,
    QStyle,
    QLineEdit,
)
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import Qt
from targets_gui import BarcodeTargetSeekerGUI  # We'll move your existing GUI here


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Barcoder Suite")
        self.setGeometry(100, 100, 600, 600)

        # Try to load icon if it exists
        icon_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "assets", "barcoder_icon.png"
        )
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))
        else:
            print(f"Icon not found at: {icon_path}")

        # Create stacked widget to handle different pages
        self.stacked_widget = QStackedWidget()
        self.setCentralWidget(self.stacked_widget)

        # Create and add pages
        self.welcome_page = self.create_welcome_page()
        self.targets_page = BarcodeTargetSeekerGUI()

        self.stacked_widget.addWidget(self.welcome_page)
        self.stacked_widget.addWidget(self.targets_page)

        # Start with welcome page
        self.stacked_widget.setCurrentWidget(self.welcome_page)

    def create_welcome_page(self):
        page = QWidget()
        layout = QVBoxLayout()

        # Welcome message
        welcome_label = QLabel("Welcome to Barcoder Suite")
        welcome_label.setAlignment(Qt.AlignCenter)
        welcome_label.setStyleSheet("font-size: 24px; margin: 20px;")

        # Try to load and display logo
        logo_path = os.path.join(
            os.path.dirname(__file__), "assets", "barcoder_logo.png"
        )
        if os.path.exists(logo_path):
            logo_label = QLabel()
            pixmap = QPixmap(logo_path)
            scaled_pixmap = pixmap.scaled(
                200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            logo_label.setPixmap(scaled_pixmap)
            logo_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(logo_label)

        # Load and display icon in the GUI
        icon_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "assets", "barcoder_icon.png"
        )
        if os.path.exists(icon_path):
            icon_label = QLabel()
            pixmap = QPixmap(icon_path)
            scaled_pixmap = pixmap.scaled(
                100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation
            )
            icon_label.setPixmap(scaled_pixmap)
            icon_label.setAlignment(Qt.AlignCenter)
            layout.addWidget(icon_label)
        else:
            print(f"Icon not found at: {icon_path}")

        layout.addWidget(welcome_label)

        # Add input field for output filename
        self.output_filename_label = QLabel("Output Filename:*")
        self.output_filename_input = QLineEdit()
        self.output_filename_input.setPlaceholderText(
            "Enter output filename (required)"
        )
        layout.addWidget(self.output_filename_label)
        layout.addWidget(self.output_filename_input)

        # Tool buttons
        targets_btn = self.create_tool_button(
            "Barcode Target Seeker",
            "Find and analyze barcode targets in your genome",
            lambda: self.stacked_widget.setCurrentWidget(self.targets_page),
        )

        layout.addWidget(targets_btn)

        # Add stretch to push everything to the top
        layout.addStretch()

        page.setLayout(layout)
        return page

    def create_tool_button(self, title, description, callback):
        button = QPushButton()
        button.clicked.connect(callback)

        # Create a container widget
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setSpacing(8)  # Increased spacing between elements
        layout.setContentsMargins(15, 15, 15, 15)  # Increased padding

        # Title with larger font and more spacing
        title_label = QLabel(title)
        title_label.setStyleSheet(
            """
            font-weight: bold;
            font-size: 18px;
            color: #212529;
            margin-bottom: 8px;
            background: transparent;
        """
        )

        # Description with better spacing
        desc_label = QLabel(description)
        desc_label.setStyleSheet(
            """
            color: #666;
            font-size: 13px;
            background: transparent;
            margin-top: 4px;
        """
        )
        desc_label.setWordWrap(True)

        layout.addWidget(title_label)
        layout.addWidget(desc_label)

        # Custom button with fixed size
        button = QPushButton()
        button.setFixedSize(400, 100)  # Set a comfortable size
        button.clicked.connect(callback)

        # Set the container as the button's content
        button.setLayout(layout)

        # Style the button
        button.setStyleSheet(
            """
            QPushButton {
                background-color: #f8f9fa;
                border: 1px solid #dee2e6;
                border-radius: 6px;
                text-align: left;
                padding: 15px;
            }
            QPushButton:hover {
                background-color: #e9ecef;
                border-color: #ced4da;
            }
        """
        )

        return button


def main():
    # Enable High DPI support
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)

    # Set application-wide stylesheet
    app.setStyleSheet(
        """
        QMainWindow {
            background-color: white;
        }
        QLabel {
            color: #212529;
        }
    """
    )

    window = MainWindow()
    window.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
