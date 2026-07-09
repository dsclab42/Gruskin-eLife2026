function saveCifti(Data, outputPath, wb_path, template,medial_mask)
Data(isnan(Data)) = 0;

num_vox_cortex = 59412;
if strcmp(outputPath(end-10),'p')
    cii_scalar = ciftiopen(template, wb_path);
    cii_scalar.cdata = Data;
    ciftisave(cii_scalar, outputPath, wb_path);
else
    cii_scalar = ciftiopen(template, wb_path);
    cii_scalar.cdata(medial_mask==1) = Data(1:num_vox_cortex);
    ciftisave(cii_scalar, outputPath, wb_path);
end
end