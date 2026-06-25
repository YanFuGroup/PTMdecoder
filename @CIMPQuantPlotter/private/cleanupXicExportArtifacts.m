function cleanupXicExportArtifacts(file_base_path)
% Delete PNG/PDF/SVG artifacts for one XIC export base path.
for ext = {'.png', '.pdf', '.svg'}
    artifact = [file_base_path, ext{1}];
    if isfile(artifact)
        delete(artifact);
    end
end
end
