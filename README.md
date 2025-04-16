# Notch_EMT

This repository contains the code for a computational modeling study investigating the interplay between Notch signaling dynamics and epigenetic regulation in melanoma metastasis. The project addresses the phenomenon where melanoma cells maintain a metastatic phenotype, potentially driven by miR-222 upregulation, even after the initiating Notch signal (triggered by cell-cell contact) is removed.

The core of the project is a dynamical systems model that integrates:
1.  **Notch Signaling Input:** Models the competition between the Notch intracellular domain (NICD) and the transcription factor MITF for binding to RBPJ.
2.  **Epigenetic Regulation of miR-222:** Simulates the histone modification states (H3K4me3 vs. H3K27me3) at the miR-222 gene locus, incorporating positive feedback loops that create bistability (memory).
3.  **Coupling:** Links Notch activity (via NICD-MITF competition) to the recruitment of histone modifying enzymes (specifically KDM5A), thus influencing the epigenetic state.

Key findings from the model analysis include:
*   The epigenetic switch mechanism enables **persistent activation** of the miR-222 locus (high H4 state) following transient Notch stimulation, providing a potential explanation for metastatic memory.
*   The system interprets **dynamic Notch signals**, acting as a **low-pass filter**. It can be robustly switched by both sustained (Dll4-like) and pulsatile (Dll1-like) signals, but higher frequency pulses require larger amplitudes or longer durations.
*   The **threshold for epigenetic switching** can be tuned by parameters within the epigenetic module, such as the feedback strength of PRC2 (H3K27me3 writer), demonstrating how intrinsic cellular factors modulate sensitivity to external signals.

Overall, this work uses computational modeling to propose and analyze a specific epigenetic mechanism underlying the persistence of Notch-induced metastatic states in melanoma, highlighting the importance of both signaling dynamics and epigenetic feedback in cell fate determination.
