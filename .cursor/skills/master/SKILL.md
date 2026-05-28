---
name: Project Master Rules
description: Python env, Django commands, English only, and the "Detective" data debugging workflow.
---

# 1. Environment & Technical Context
- **Python Path**: `/Users/nht435/miniconda3/envs/gpcrdb/bin/python`
  - Always use this interpreter.
- **Project Type**: Django (GPCRdb).
- **Start Command**: `python manage.py runserver 0.0.0.0:8000`

# 2. Language & Communication
- **Primary Language**: Strictly **English** for code, logic, and explanations.
- **Handling Chinese Input**:
  - If I ask in **Chinese**: You must answer in **English**, but provide a **brief Chinese summary** at the end.
- **Code Comments**: English only.

# 3. Tone & Style
- Be professional, concise, and direct.

# 4. Data Detective Protocol (Debugging Strategy)
When investigating data discrepancies or logic bugs, strictly follow this workflow:

## A. Execution Method
- **Never** put debug print statements into `views.py` or `models.py`.
- **Always** create a standalone script (e.g., `debug_investigation.py`).
- **Run Command**: `python manage.py shell < debug_investigation.py`

## B. Script Structure (The "Report" Style)
The script must produce a readable "Detective Report" following this template:
1.  **Imports**: Include necessary models and tools (`from django.db.models import F, Q, Count`).
2.  **Header**: Use `print("="*60)` and a descriptive title (e.g., "Detective Report: Missing Models").
3.  **Step-by-Step Logic**:
    - **Step 1: The Baseline**: Query the data based on current logic. Print the count.
    - **Step 2: The Hypothesis**: Query based on an alternative/older logic or wider scope. Print the count.
    - **Step 3: The Gap**: Calculate and print the difference (`count_hypothesis - count_baseline`).
4.  **Output Formatting**:
    - Use f-strings for clarity: `print(f"1. Current Logic count: {count}")`.
    - If lists are long, only print the first 5 IDs: `ids[:5]`.