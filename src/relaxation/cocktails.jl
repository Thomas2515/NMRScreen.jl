function prepcocktails(cocktails)
    smiles = []
    for cocktail in cocktails
        for peak in cocktail.peaks
            push!(smiles, peak.smiles)
        end
    end
    embeddings = get_embeddings(smiles)

    # input: cocktails where peaks are BasicPeaks
    # output: cocktails where peaks are ScreeningPeaks, also contains UMAP coordinates
    for cocktail in cocktails
        for (i, peak) in enumerate(cocktail.peaks)
            newpeak = ScreeningPeak(peak, embeddings[1, i], embeddings[2, i])
            cocktail.peaks[i] = newpeak
        end
        sort!(cocktail.peaks, by=sp -> sp.reference_shift, rev=true)
    end

    return cocktails
end

function get_embeddings(smiles)
    fingerprint_matrix = []
    for i in smiles
        mol = get_mol(i)
        fingerprint = get_morgan_fp(mol)
        fingerprint_split = split(fingerprint, "")
        fingerprint_parsed = []
        for bits in fingerprint_split
            push!(fingerprint_parsed, parse.(Float64, bits))
        end
        push!(fingerprint_matrix, fingerprint_parsed)
    end
    fingerprint_matrix = hcat(fingerprint_matrix...)
    return umap(fingerprint_matrix, 2; n_neighbors=20, min_dist=0.1, metric=Jaccard())
end
