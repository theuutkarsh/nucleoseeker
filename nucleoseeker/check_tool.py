import subprocess
import sys

def check_tool(tool_name):
    """Check if a tool is installed and available in the system PATH."""
    try:
        subprocess.run([tool_name, '-h'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        sys.exit(f"Error: {tool_name} is not installed. Please install it to proceed.")
