import { defineConfig } from 'astro/config';
import { unified } from '@astrojs/markdown-remark';
import starlight from '@astrojs/starlight';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';
import navigation from './navigation.json' with { type: 'json' };

function mapSidebarItems(items) {
  return items.map((item) => {
    if (item.type === 'page') {
      return { slug: item.slug, ...(item.label ? { label: item.label } : {}) };
    }

    if (item.type === 'group') {
      return {
        label: item.label,
        ...(item.collapsed === true ? { collapsed: true } : {}),
        items: mapSidebarItems(item.items),
      };
    }

    return {
      label: item.label,
      link: item.link,
      ...(item.newTab
        ? { attrs: { target: '_blank', rel: 'noreferrer' } }
        : {}),
    };
  });
}

// A section with `autogenerate` lists generated pages, such as the manuscript chapters.
const sidebar = navigation.sections.map((section) => ({
  label: section.label,
  ...(section.autogenerate
    ? { items: [{ autogenerate: { directory: section.autogenerate } }] }
    : { items: mapSidebarItems(section.items) }),
}));

export default defineConfig({
  site: 'https://nkzono99.github.io',
  base: '/BEACH',
  outDir: '../build/starlight-site',
  build: {
    format: 'file',
  },
  markdown: {
    processor: unified({
      remarkPlugins: [remarkMath],
      rehypePlugins: [[rehypeKatex, { strict: false, throwOnError: true }]],
    }),
  },
  integrations: [
    starlight({
      title: 'BEACH',
      defaultLocale: 'root',
      // Japanese is the canonical language while the docs are restructured (ADR 0003).
      locales: {
        root: {
          label: '日本語',
          lang: 'ja',
        },
      },
      customCss: ['./src/styles/custom.css'],
      lastUpdated: true,
      components: {
        Head: './src/components/Head.astro',
      },
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/Nkzono99/BEACH',
        },
      ],
      sidebar,
    }),
  ],
});
