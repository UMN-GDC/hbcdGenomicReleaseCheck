import os
import xml.etree.ElementTree as ET
import requests

# 1. Configuration
GIST_ID = "YOUR_COPIED_GIST_ID"
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")  # Kept secure on your local machine
XML_FILE = "reports.xml"

if not GITHUB_TOKEN:
    raise ValueError("Please set the GITHUB_TOKEN environment variable.")

# 2. Parse the local JUnit XML file
try:
    tree = ET.parse(XML_FILE)
    root = tree.getroot()
    testsuite = root if root.tag == "testsuite" else root.find("testsuite")

    total = int(testsuite.attrib.get("tests", 0))
    failures = int(testsuite.attrib.get("failures", 0))
    errors = int(testsuite.attrib.get("errors", 0))
    passed = total - failures - errors

    if failures == 0 and errors == 0:
        status_msg = f"{passed} passed"
        color = "green"
    else:
        status_msg = f"{failures} failed, {passed} passed"
        color = "red"

except Exception as e:
    # If the file is missing or broken, default to an error status
    status_msg = "error running tests"
    color = "inactive"

# 3. Construct the Shields.io JSON payload
payload = {
    "files": {
        "test-status.json": {
            "content": f'{{"schemaVersion": 1, "label": "tests", "message": "{status_msg}", "color": "{color}"}}'
        }
    }
}

# 4. Push the payload up to your GitHub Gist
headers = {
    "Authorization": f"Bearer {GITHUB_TOKEN}",
    "Accept": "application/vnd.github+json",
}
response = requests.patch(
    f"https://api.github.com/gists/{GIST_ID}", json=payload, headers=headers
)

if response.status_code == 200:
    print("Successfully updated GitHub landing page badge!")
else:
    print(f"Failed to update Gist: {response.status_code}, {response.text}")
