# SMLapSVM
Multicategory Laplacian support vector machines

```SMLapSVM``` is an R package. ```SMLapSVM``` provides functions for fitting the reinforced multicategory Laplacian support vector machine (RMLapSVM) and the reinforced angle-based multicategory Laplacian support vector machine (RAMLapSVM).
Furthermore, ```SMLapSVM``` also provides functions to fit the structured RMLapSVM (SRMLapSVM) and the structured RAMLapSVM (SRAMLapSVM). 

## 1. INSTALLATION

(1) From GitHub
```{r}
> library(devtools)
> install_github("bbeomjin/SMLapSVM")
```

## 2. USAGE NOTES

(1) Description of R functions in ```SMLapSVM```

- Descriptions of arguments in the functions in ```SMLapSVM``` can be obtained by help() or ? in R prompt, and documentation of ```SMLapSVM```.   

(2) List of R functions in ```SMLapSVM``` package

- ```mlapsvm``` : A function for fitting the RMLapSVM and RAMLapSVM.

- ```smlapsvm``` : A function for learning the SRMLapSVM and SRAMLapSVM.

- ```cv.mlapsvm``` : 
