\* Problem: MFPsub *\

Maximize
 obj: - 276 y2

Subject To
 LB0: + L0 = 0
 LB1: + L1 - L0 = 0
 LB2: - 12 y2 + L2 - L1 = 0
 LB3: + L3 - L2 = 0

Bounds
 0 <= L0 <= 300
 0 <= L1 <= 300
 0 <= L2 <= 300
 0 <= L3 <= 300
 0 <= y1 <= 1
 0 <= y2 <= 1
 0 <= y3 <= 1

End
