# Workstationinstruktioner

##### Til opsætning af sekventering

Når biblioteket er færdiglavet, og flowcellen skal til at pakkes ud, kan denne opskrift bruges til at starte selve sekventeringen. Hvis der opstår problemer, eller computeren opfører sig sært, kontakt da Carl på tlf. 5117 5719.



1. <u>Inden biblioteket loades på flowcellen udføres et **flow-cellecheck**.</u>
   1. Log in på en ledig ontseq-workstation; login-oplysningerne står på en seddel som hænger på skærmen.
   2. Åbn **MinKNOW**, ved at klikke på ikonet i venstre side af skærmen, eller tryk på windows-tasten og søg efter minknow.
      1. Hvis MinKNOW beder om et login, så brug det som står på den seddel der hænger på skærmen.
   3. Inde i MinKNOW åbnes nu **Start**, og herefter vælges **Flow cell check**. Denne process tager noget tid, og skulle gerne rapportere antallet af aktive porer på cellen. 
2. <u>Efter et successfuldt flow-cellecheck, skal flow-cellen nu loades.</u>
   1.  Følg instruktionerne omkring loading i "den anden" manual.
   2. Indsæt flowcellen i en Minion-sekvenator.
3. <u>Når flow-cellen er loadet med bibliotek skal sekventeringen sættes i gang.</u>
   1. Tilkobl sekvenatoren med et micro-USB3 til USB3-A kabel.
      1. Verificer, at USB3-A enden er isat i et "SuperSpeed" (SS) stik i computeren.
   2. Åbn **MinKNOW** igen, som vi tidligere brugte til at lave flow-cellecheck med.
   3. Vælg **Start** og derefter **Start sequencing**.
      1. Side: **1. Positions**
         1. På side **1. Positions** skrives `YYYY-MM-DD` (fx `2021-03-15`) i **Experiment name**-feltet. Da det tager mindst et døgn at køre et en sekventering, er der ikke fare for duplikater, og det er derfor tilstrækkeligt at bruge datoen alene, som batch-identifikation. 
         2. Verificer at der i **Positions**-tabellen står "FLO-MIN106" under **Flow cell type**
         3. Feltet **Sample ID** skal forblive tomt.
         4. Tryk på **Load saved settings** og vælg "ontseqX SARS-CoV-2" hvor X er maskinens ID (A eller Z). Maskinens navn kan ses på det klistermærke der står på skærmen. Tryk på **Load**.
            1. Der kommer en pop-up som spørger om du ønsker at loade kit "SQK-LSK109", vælg **Yes**.