## Instruktioner til opsætning af sekventering på sekventerings-workstations

Når biblioteket er færdiglavet, og flowcellen skal til at pakkes ud, kan denne opskrift bruges til at starte selve sekventeringen. Hvis der opstår problemer, eller computeren opfører sig sært, kontakt da Carl på tlf. 5117 5719 eller carkob@rm.dk



1. <u>Inden biblioteket loades på flowcellen udføres et **flow-cellecheck**.</u>
   1. Log ind på en ledig ontseq-workstation; login-oplysningerne står på en seddel som hænger på skærmen.
   2. Åbn **MinKNOW**, ved at klikke på ikonet i venstre side af skærmen, eller tryk på windows-tasten og søg efter minknow.
      1. Hvis MinKNOW beder om et login, så brug det som står på den seddel der hænger på skærmen.
   3. Inde i MinKNOW åbnes nu **Start**, og herefter vælges **Flow cell check**. Denne process tager noget tid, og skulle gerne rapportere antallet af aktive porer på cellen. 
2. <u>Efter et successfuldt flow-cellecheck, skal flow-cellen nu loades.</u>
   1.  Følg instruktionerne omkring loading i "den anden" manual.
   2. Indsæt flowcellen i en Minion-sekvenator.
3. <u>Når flow-cellen er loadet med bibliotek skal sekventeringen sættes i gang.</u>
   1. Tilkobl sekvenatoren med det kabel der stikker ud af skærmen.
   2. Åbn **MinKNOW** igen, som vi tidligere brugte til at lave flow-cellecheck med.
   3. Vælg **Start** og derefter **Start sequencing**.
      1. På side **1. Positions** skrives `YYYY-MM-DD` (fx `2021-03-15`) i **Experiment name**-feltet.
      3. Feltet **Sample ID** skal forblive tomt.
      5. Tryk på **Load saved settings** og vælg "ontseqX SARS-CoV-2" hvor X er maskinens ID ("a" eller "z"). Maskinens navn kan ses på den seddel der sidder på skærmen. Tryk på **Load**.
         1. Der kommer en pop-up som spørger om du ønsker at loade kit "SQK-LSK109", vælg **Yes**.
      6. Vælg **Skip to end >>** og tjek at alle instillingerne er som på den printede seddel som hænger på skærmen (TODO: ikke printet i skrivende stund.)
      7. Tryk til sidst på **Start**-knappen nede i højre hjørne.
      8. Gå nu ind i **Experiments** (venstre menu) og hold øje med **Run state** for den startede kørsel.
      9. Du vil se meddelelserne: "_Heating to N degrees_" efterfulgt af "_Performing MUX scan_".
      8. Gå ind i **System messages (venstre menu)** og tjek at der ikke er nogle advarsler eller fejlmeddelelser.
         1. Hvis der kommer en advarselsmeddelelse om at der ikke er så meget plads tilbage på harddisken *"Disk usage alert"*, kan du tjekke at der er mindst 500 GB ledig plads. Dette gøres ved at åbne **System Monitor** ikon under terminal-ikonet i venstre side af skærmen, og gå ind under **File Systems**-fanen.
      9. Når **Run state** viser "_Active_", kan du klikke ind på kørslens undersider. Her kan du se generelle statistikker omkring kørslen, og få en ide om hvor mange porer der er aktive.
   4. Hvis du er tilfreds med kørslens statistikker, kan du starte **Pappenheim**-pipelinen. Pappenheim starter **Rampart** som er et sekventerings-overvågningsværktøj der hjælper med at danne et overblik over hvordan sekventeringen skrider frem. Pappenheim mapper efter endt sekventering reads til et SARS-CoV-2 referencegenom, således at vi på længere sigt kan kalde varianter. Pappenheim startes i terminalen og er udelukkende tekstbaseret. 
      1. Før du starter pappenheim, skal du have et *sample sheet* i xlsx-format klar. Dette *sample sheet* skal som minimum indeholde to kolonner: **Barcode** og **Sample ID**. Det er ligemeget hvilken rækkefølge kolonnerne er i, og det er også ligemeget hvilke andre kolonner der er tilstede i dit *sample sheet*. Det eneste der er vigtigt er, at pappenheim får mulighed for at regne ud hvilke barcodes der passer sammen med hvilke prøver. Sekventerings-computerene har adgang til KMAs delte drev, hvorfra du kan tilgå dit samplesheet. Ellers kan du overføre dit *sample sheet* til workstationen enten ved at bruge et USB-stik eller ved at sende dig selv en mail og hente den via webmail i computerens browser. Når du har dit *sample sheet* klar på computeren, er du klar til at fortsætte til næste trin:
      2. Åbn et **Terminal**-vindue ved at klikke på **Terminal**-ikonet i venstre side af skærmen. Ikonet skulle gerne ligge lige under MinKNOW-ikonet. Hvis ikke du kan finde dette ikon, så tryk på windows-tasten (på tastaturet) og søg efter *terminal*. I terminal-prompten skriver du **start_pappenheim** efterfulgt af et mellemrum. Dernæst skal du angive stien til dit samplesheet. Igen efterfulgt af et mellemrum skal du angive kørsels-stien ("rundir") hvor minknow er blevet bedt om at outputte sekvensdata. Den oplagte måde at angive samplesheet og kørsels-sti er ved at finde samplesheetet i stifinderen og trækker "drag-and-drop" over i terminal vinduet. Det samme gøres for kørsels-stien.
         *. Her ses et generelt eksempel:
         ![alt text](https://github.com/KMA-Aarhus/pappenheim/blob/main/documentation/generelt%20eksempel.png)
         
         *Genvej til at finde samplesheet: Åbn stifinderen, og vælg genvejen i venstre spalte som hedder "NanoPore". Heri har du måske allerede lagt dit samplesheet*
         *Genvej til at finde kørsels-stien: Åbn stifinderen, og vælg genvejen i venstre spalte som hedder "sc2_sequencing". MinKNOW skulle gerne være indstillet til automatisk at skrive rådata hertil.*
         
         Når du har udfyldt kommandoprompten, kan du trykke 'enter' på tastaturet, og pappenheim-pipelinen skulle gerne sørge for resten.
         



Hvis der opstår problemer er du meget velkommet til at kontakte Carl på carkob@rm.dk eller +45 51175719
