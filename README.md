# SpatialDE_app
Identify spatially variable genes based on the input data of spatial cells positions, and cells genes expression for Merfish, exseq, smFish.

### How to run this:
1. Prepare inputs: 
    - "gene by cell" count matrix,
    -  metadata with spatial coordinate columns named as "center_x" and "center_y" in the same cell order.
    
2. Run the command in terminal:
    ```python
    python cli.py run -c [COUNT MATRIX PATH] -m [METADATA PATH]
    ```