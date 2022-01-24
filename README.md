# Exseq Analysis
Identify spatially variable genes based on the input data of spatial cells positions, and cells genes expression for exseq.

### Analysis Pipeline and Results
- `spatialDE_analysis.py` : results are in `220116_SpatialDE_results.zip`
- `hotspot_analysis.py` : results are in `220123_Hotspot_results.zip`


### How to run SpatialDE:
1. Prepare inputs: 
    - "gene by cell" count matrix,
    -  metadata with spatial coordinate columns named as "center_x" and "center_y" in the same cell order.
    - optional: batch correct needed or not needed, use --batch-correct.
    
2. Run the command in terminal:
    ```python
    python cli.py run -c [COUNT MATRIX PATH] -m [METADATA PATH]
    ```
