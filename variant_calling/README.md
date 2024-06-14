# Variant Calling Pipeline/Workflow

*Scripts in /original_scripts were created by Natalie Gonzales. We have changed various things across the years to optimize for our respective projects.*

For the most part, all three workflows are the same. The only difference is the preprocessing step. The Legacy workflow uses the original preprocessing script, while the Alternative workflow uses a modified version of the original preprocessing script. The Troubleshooting/diagnostic workflow uses the original preprocessing script, but adds additional steps to diagnose and fix issues.

The reason to "split" into these three workflows is to facilitate troubleshooting of failed .bam files. Some of the errors identified are lack of proper read group assignment, unpaired mates and missing mate information. The troubleshooting/diagnostic workflow adds steps to address these issues while the Alternative workflows aims to address the issues from the start. 

## Ferris Legacy workflow:

```mermaid
graph TD
    A[Step 1] --> B[Step 2_Legacy] --> C[Step 3] --> D[Step 4] 
```
### The scripts for the Legacy workflow are:
```
variant_calling/
    ├── Step1_readQC.sh
    ├── Step2_Legacy_preprocessing_simple.sh
    ├── Step3_HaplotypeCaller_simple.sh
    ├── Step4_JointGenotyping_simple.sh
```

## Troubleshooting/diagnostic workflow:

```mermaid
graph TD
    A[Step 1] --> B[Step 2] --> T
        subgraph Diagnostic
        T(Step 2_1) --> X(Step 2_2) --> Y(Step 2_3)
    end
 B .-> C
    Y --> C
    C[Step 3] --> D[Step 4]    
```
### The scripts for the Troubleshooting/diagnostic workflow are:
```
variant_calling/
    ├── Step1_readQC.sh
    ├── Step2_1_diagnostics.sh
    ├── Step2_2_ReassigningRG.sh
    ├── Step2_3_FixMateInformation.sh
    ├── Step2_Legacy_preprocessing_simple.sh
    ├── Step3_HaplotypeCaller_simple.sh
    ├── Step4_JointGenotyping_simple.sh
```

## Alternative workflow:

```mermaid
graph TD
    A[Step 1] --> B[Step 2_Alter] --> C[Step 3] --> D[Step 4]
```
### The scripts for the Alternative workflow are:
```
variant_calling/
    ├── Step1_readQC.sh
    ├── Step2_Alter_preprocessing.sh
    ├── Step3_HaplotypeCaller_simple.sh
    ├── Step4_JointGenotyping_simple.sh
```