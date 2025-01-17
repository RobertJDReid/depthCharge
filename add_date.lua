-- add_date.lua
function Meta(meta)
  if meta.date == nil then
    meta.date = os.date("%B %e, %Y")
    return meta
  end
end
