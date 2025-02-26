ChemicalizeMarvinJs.createEditor("#marvin-editor").then(function (marvin) {
    marvin.setDisplaySettings({
        "toolbars": "education",
    });
    function transferStructure() {
        marvin.exportStructure('cxsmiles', {'extra': 'f'}).then(function (smiles) {
            $("#smiles").val(smiles)
        })
    }
    document.getElementById("transfer-structure").addEventListener("click", transferStructure);
})