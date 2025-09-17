# Core Code Changes for MaBoSS Mutation Tracking in PhysiCell

This document summarizes the key changes and additions made to the core PhysiCell and PhysiBoSS codebase to enable random node mutation, mutation tracking, and lineage tracing.  
**Files referenced:**  
- `addons/PhysiBoSS/src/maboss_intracellular.cpp`
- `addons/PhysiBoSS/src/maboss_intracellular.h`
- `addons/PhysiBoSS/src/maboss_network.cpp`
- `addons/PhysiBoSS/src/maboss_network.h`
- `modules/PhysiCell_MultiCellDS.cpp`
- `modules/PhysiCell_MultiCellDS.h`
- `core/PhysiCell_cell.cpp`
- `core/PhysiCell_cell.h`
- `custom_modules/custom.cpp` (for the mutation function)

---

## 1. **Random Node Mutation and Mutation Logging**

**File:** `custom_modules/custom.cpp`  
**Function:** `phase_exit_mutation_function`

- Added logic to:
  - Select a random node from the MaBoSS network for a cell.
  - Flip its boolean state.
  - Record the mutation as a string (e.g., `NodeName_1`) in the cell's `mutations` vector.
- Ensured the function checks for the presence of a MaBoSS model before acting.

**Key code:**
```cpp
if (pCell->phenotype.intracellular &&
    pCell->phenotype.intracellular->intracellular_type == "maboss")
{
    MaBoSSIntracellular* maboss_model = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular);
    std::vector<std::string> node_names = maboss_model->maboss.get_all_node_names();
    // ... select random node, flip, and log mutation ...
}
```

---

## 2. **Mutation Vector Output Formatting**

**File:** `modules/PhysiCell_MultiCellDS.cpp`

- Changed the output of the `mutations` vector in the `.mut.csv` file to use `"."` as a separator instead of spaces.
- Ensured the output is consistent for downstream analysis.

**Key code:**
```cpp
std::string buffer = "";
for( int j=0 ; j < pCell->custom_data.mutations.size(); j++ )
{
    buffer += pCell->custom_data.mutations[j];
    if (j != pCell->custom_data.mutations.size() - 1)
        buffer += ".";
}
file_mut << buffer << std::endl;
```

---

## 3. **Cell Lineage Tracking**

**Files:**  
- `core/PhysiCell_cell.h`
- `core/PhysiCell_cell.cpp`

- Added two new attributes to the `Cell` class:
  - `int generation` (starts at 0 for initial cells, increments on division)
  - `int parent_ID` (tracks the parent cell's ID, -1 for initial cells)
- Set these attributes during cell creation and division.

**Key code in `Cell` class:**
```cpp
int generation = 0;   // Track cell generation (0 for initial cells)
int parent_ID = -1;   // Track parent cell ID (-1 for initial cells)
```
**Set in `Cell::divide()`:**
```cpp
child->generation = this->generation + 1;
child->parent_ID = this->ID;
```

---

## 4. **Mutation CSV Output Enhancements**

**File:** `modules/PhysiCell_MultiCellDS.cpp`

- Extended the `.mut.csv` output to include:
  - Simulation time
  - Cell ID
  - Cell type
  - Parent ID
  - Generation
  - Mutations (dot-separated)

**CSV header:**
```
time,ID,type,parent_ID,generation,mutations
```

**CSV row example:**
```
120.0,42,1,17,2,A_1.B_0
```

---

## 5. **MaBoSS API Usage**

**Files:**  
- `addons/PhysiBoSS/src/maboss_intracellular.cpp`
- `addons/PhysiBoSS/src/maboss_intracellular.h`
- `addons/PhysiBoSS/src/maboss_network.cpp`
- `addons/PhysiBoSS/src/maboss_network.h`

- Used the public MaBoSS API to:
  - Retrieve all node names: `get_all_node_names()`
  - Get/set node values: `get_node_value(name)`, `set_node_value(name, value)`
- No changes to these files were required for the above API, but ensure your version includes these methods.

---

## 6. **Include Statements**

- In any file using `MaBoSSIntracellular`, ensure the following include is present:
  ```cpp
  #include "../addons/PhysiBoSS/src/maboss_intracellular.h"
  ```

---

## 7. **Summary Table**

| Feature/Change                | File(s) Affected                        | Description                                 |
|-------------------------------|-----------------------------------------|---------------------------------------------|
| Random node mutation & logging| `custom.cpp`                            | Flip random node, log mutation              |
| Dot-separated mutation output | `PhysiCell_MultiCellDS.cpp`             | Output mutations as `A_1.B_0`               |
| Lineage tracking              | `PhysiCell_cell.h/cpp`                  | Add `generation`, `parent_ID` to `Cell`     |
| Extended mutation CSV         | `PhysiCell_MultiCellDS.cpp`             | Add time, parent_ID, generation to CSV      |
| MaBoSS API usage              | `maboss_intracellular.*`, `maboss_network.*` | Use public API for node access         |
| Required includes             | Any using MaBoSS types                  | Add correct `#include`                      |

---

**Keep this file up to date as you make further changes!**
